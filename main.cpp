#include <iostream>
#include <fstream>
#include <string>
#include <gdal.h>
#include <gdal_priv.h>
#include <map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <limits>
#include <list>
#include <array>
#include <chrono>

struct meta{
        std::string file;
        int rows;
        int cols;
        GDALDataType gdal_type;
        int byte_size;
        meta(std::string _file) : file(_file) {};
        meta() {};
        friend std::ostream& operator << (std::ostream& os,const meta &m);
        // meta& operator=(meta &&copy){
        //     // file = std::move(copy.file);
        //     return *this;
        // }

    };

    std::ostream &operator<<(std::ostream &os,const meta& m){
            os << "file: " << m.file << " / ncols: "  <<  m.cols << " / nrows: " << m.rows << " / type: " << m.gdal_type << std::endl;
           return os;
}

enum class Direction : uint {
    No_outlet = 0,
    Northeast = 1,
    East = 2,
    Southeast = 4,
    South = 8,
    Southwest = 16,
    West = 32,
    Northwest = 64,
    North = 128
};    

struct Orientation {
    enum class quadrant  {
        I = 0,
        II = 1,
        III = 2,
        IV = 3,
        UP = 4,
        RIGHT = 5,
        DOWN = 6,
        LEFT = 7 
    };

    Orientation (int dx,int dy){
    if      (dx == 0 && dy  > 0) current_orientation = quadrant::UP;
    else if (dx == 0 && dy  < 0) current_orientation = quadrant::DOWN;
    else if (dx  > 0 && dy == 0) current_orientation = quadrant::RIGHT;
    else if (dx  < 0 && dy == 0) current_orientation = quadrant::LEFT;
    else if (dx  > 0 && dy  > 0) current_orientation = quadrant::I;
    else if (dx  > 0 && dy  < 0) current_orientation = quadrant::II;
    else if (dx  < 0 && dy  < 0) current_orientation = quadrant::III;
    else if (dx  < 0 && dy  > 0) current_orientation = quadrant::IV;
    }

    quadrant current_orientation;

    void transform(std::list<std::pair<int,int>>& coord_list){
        switch (current_orientation){
            case quadrant::DOWN: 
                for (auto&i: coord_list) i.second *= -1;
                break;
            case quadrant::RIGHT:
                for (auto&i: coord_list) i.first *= -1;
                break;
            case quadrant::II:
                for (auto&i: coord_list) i.second *= -1;
                break;
            case quadrant::III:
                for (auto&i: coord_list) {
                    i.second *= -1; 
                    i.first *= -1;
                }
                break;
            case quadrant::IV:
                for (auto&i: coord_list) i.first *= -1;
                break;
            default:
                break;
        }; 
     }

};

struct pix{
    float elevation;
    Direction dir; 
    size_t flowacc;
    size_t idx;
    pix* downstream;
    std::vector<pix*>upstream;
    pix* max_upstream;
};

enum Condition {Upstream,River,Dam,Dam_curtain};     

class Up_condition{
    public:
        virtual bool check(const pix* up, const pix* down) = 0;
        virtual ~Up_condition() = default; 
        static bool follow_river;
        static pix* outlet;
};

bool Up_condition::follow_river = false;
pix* Up_condition::outlet = nullptr;

class Up_condition_upstream: public Up_condition {
    public:
        bool check(const pix* up, const pix* down){
            return (up->downstream == down);
        }
};

class Up_condition_dam: public Up_condition {
    public:
        bool check(const pix* up, const pix* down){
            return (up->downstream == down && (up->elevation - outlet->elevation) <= 4.0);
        }
}; 

class Up_condition_dam_curtain : public Up_condition{
    public :
    bool check (const pix* up, const pix* down){
        return true;
    }
};

class Up_condition_river: public Up_condition {
    public:
        bool check(const pix* up, const pix* down){
            num_iterations++;
            if (num_iterations < max_num_upstream_pixels) return true; 
            else {
                num_iterations = 0;
                return false;
            } 
        }
        Up_condition_river() {follow_river=true;};
        uint max_num_upstream_pixels = 100;
        uint num_iterations = 0;
};

class Creator {
    public:
        virtual std::unique_ptr<Up_condition> create(Condition cond){
            if(Condition::Upstream == cond) return std::make_unique<Up_condition_upstream>();
            if(Condition::Dam == cond) return std::make_unique<Up_condition_dam>();
            if(Condition::River == cond) return std::make_unique<Up_condition_river>();
            if(Condition::Dam_curtain == cond) return std::make_unique<Up_condition_dam_curtain>();
            return nullptr;
        }
        virtual ~Creator() = default;        
};

class Data{
    public:

        // friend std::ostream& operator << (std::ostream& os,const Data &m);
        Data(){};
        // struct idx_comp {
        //     bool operator () (const Data& d, const size_t idx) const {return d.idx < idx;}
        //     bool operator () (const size_t& idx, const Data& d) const {return idx < d.idx;}
            
        // };

        static void set_band_size(int _nrows, int _ncols){
            nrows = _nrows;
            ncols = _ncols;
            numel = (size_t) nrows * (size_t) ncols;
        }

        static void set_number_valid_elements(size_t num){
            num_valid_elements =  num;
        }

        static std::pair<int,int> get_xy(size_t& index){
            std::pair<int,int> result;
            result.first = index % ncols;
            result.second = index / ncols;
            return result;
        }

        static void set_downstream_idx_for_all_data(){
            for (size_t i = 0; i < num_valid_elements;++i){
                size_t downstream_idx = get_index_in_direction(all_data[i].idx,all_data[i].dir);
                if (downstream_idx == no_value){
                    all_data[i].downstream=nullptr;
                }
                else{
                    all_data[i].downstream = binary_search_pix(downstream_idx);
                }

            }
        }

        static void set_upstream_with_largest_flowacc(){
            for (auto i = 0; i < num_valid_elements; ++i){
                pix* pt;
                size_t max_accumulation = 0;
                if (! all_data[i].upstream.empty()) {
                    for (std::vector<pix*>::iterator j = all_data[i].upstream.begin(); j !=all_data[i].upstream.end(); ++j) {
                        if ((*j)->flowacc > max_accumulation) {
                            pt = *j;
                            max_accumulation = (*j)->flowacc;
                        }
                    }
                } else pt = nullptr;
                
                all_data[i].max_upstream = pt;
            }
        }

        static size_t get_index_in_direction(size_t idx, Direction dir){
            auto coords =  get_xy(idx);
            switch (dir){
                case Direction::No_outlet :
                    return no_value;
                    // break;
                case Direction::North:
                    coords.second -= 1;
                    break;
                case Direction::South:
                    coords.second += 1;
                    break;
                case Direction::East:
                    coords.first += 1;
                    break;
                case Direction::West:
                    coords.first -= 1;
                    break;
                case Direction::Northwest:
                    coords.second -= 1;
                    coords.first -= 1;
                    break;
                case Direction::Northeast:
                    coords.second -= 1;
                    coords.first += 1; 
                    break;
                case Direction::Southwest:
                    coords.second +=1;
                    coords.first -=1;
                    break;
                case Direction::Southeast:
                    coords.second += 1;
                    coords.first += 1;
                    break;    
            }

            if (coords.first < 0 || coords.first >= ncols || coords.second < 0 || coords.second >= nrows){
                return no_value;
            }
            else {
                return (coords.second * ncols + coords.first);
            }
        }

        static const size_t no_value = std::numeric_limits<size_t>::max();

        void get_basin(pix* outlet){
            basin.clear();
            basin.push_back(outlet);
            cond->outlet = outlet;
            get_all_upstream(outlet);
            std::cout << "Number of upstream pixels for this pixel: " << basin.size() << std::endl;
        }

        void get_all_upstream(pix* data){ 
            if (!cond->follow_river){
                for (auto& i : data->upstream){
                    if (cond->check(i,data)){
                        basin.push_back(i);
                        // if (cond->update_outlet) data = i;
                        get_all_upstream(i);
                    }                
                } 
            }
            else {
                basin.push_back(data->max_upstream);
                auto coord = get_xy(data->idx);
                // std::cout << data->idx << " " << coord.first << " " << coord.second << " " <<  data->flowacc << std::endl;
                if (data->max_upstream != nullptr && cond->check(nullptr,nullptr)) get_all_upstream(data->max_upstream);                 
            }
                      
        }

        std::vector<pix*> get_curtain(pix* outlet, pix* river_upstream, uint max_width=20, uint max_height = 2){
            auto start_coord = get_xy(outlet->idx);
            auto end_coord = get_xy(river_upstream->idx);
            int dx = end_coord.first - start_coord.first;
            int dy = end_coord.second - start_coord.second;
            double slope = std::numeric_limits<double>::max();
            if (dx != 0) slope = (double) dy/ (double) dx;
            //extending river
            int current_distance = dx * dx + dy * dy;
            int multiplier = 2;
            if ( !(current_distance > (multiplier * (max_width * max_width) ) ) ){
                dx *= multiplier;
                dy *= multiplier;
            }
            //Getting perpendicular pixels
            auto clockwise_perpendicular = start_coord;
            clockwise_perpendicular.first += dy;
            clockwise_perpendicular.second -= dx;
            auto counterclockwise_perpendicular = start_coord;
            counterclockwise_perpendicular.first -= dy;
            counterclockwise_perpendicular.second += dx;
            
            //Tracing perpendicular pixels
            std::cout << "Clockwise" << std::endl;
            auto clockwise_coordinates = bresenham(start_coord,clockwise_perpendicular);
            for (auto&i : clockwise_coordinates) std::cout << "(" << i.first << "," << i.second << ")" << std::endl;
            std::cout << "Counterclockwise" << std::endl;
            auto counterclockwise_coordinates = bresenham(start_coord,counterclockwise_perpendicular);
            for (auto&i : counterclockwise_coordinates) std::cout << "(" << i.first << "," << i.second << ")" << std::endl;

            std::list<pix*> clockwise_pixels, counterclockwise_pixels;
            pix* center = binary_search_pix(get_idx(clockwise_coordinates.begin()->first,clockwise_coordinates.begin()->second)); 
            clockwise_coordinates.pop_front();
            counterclockwise_coordinates.pop_front();
            for (auto& p : clockwise_coordinates) clockwise_pixels.push_back(binary_search_pix(get_idx(p.first,p.second))); 
            for (auto& p : counterclockwise_coordinates) counterclockwise_pixels.push_back(binary_search_pix(get_idx(p.first,p.second))); 
            
            float valid_height = center->elevation + (float) max_height; 
            uint squared_dist = max_width * max_width;
            std::vector<pix*> result;
            result.push_back(center);
            for(auto&p : clockwise_pixels){
                if (p->elevation <= valid_height ){
                    result.push_back(p);
                }
            }

            for(auto&p : counterclockwise_pixels){
                if (p->elevation <= valid_height ){
                    result.push_back(p);
                }
            }

            return result;

            //saving pixel coordinates to file
            // std::ofstream out;
            // out.open("dam.txt");
            // for (auto it = clockwise_pixels.rbegin(); it != clockwise_pixels.rend(); ++it) out << it->first << "," << it->second << "\n";
            // auto it = counterclockwise_pixels.begin()++;
            // for(;it!=counterclockwise_pixels.end();++it) out << it->first << "," << it->second << "\n";
            // out.close();

            //Getting the pixels that compose the curtain
            
            

        }

        static void set_all_upstream_pixels(){
            for (auto i = 0; i < num_valid_elements; ++i){
                if (all_data[i].downstream != nullptr){
                    all_data[i].downstream->upstream.push_back(&all_data[i]);
                } 
                // pixels[i].set_as_upstream_pixel();
            }
        }

        
        static size_t get_idx(int x, int y) {return  (size_t) y * (size_t) ncols + (size_t) x;}


        static pix* binary_search_pix(size_t idx_to_find){
            size_t l = 0;
            size_t r = num_valid_elements;

            if (idx_to_find == Data::no_value){
                return nullptr;
            }
            else {
                while  (l != r ){
                    size_t m = l + ( r - l )/2;
                    if (idx_to_find < all_data[m].idx){
                        r = m;
                    }
                    
                    else if (all_data[m].idx < idx_to_find) {
                        l = m + 1;
                    }
                    else {
                        return &all_data[m];
                    }               
                }
                return nullptr;
            }
        }
    

        static void set_all_data(pix* _all_data){
            all_data = _all_data;
        }
    
        std::vector<pix*> basin;
        void paint_pixels(uint8_t* save_data){
            for (auto &p: basin) save_data[p->idx] = 1;
        }


    static std::list<std::pair<int,int>>  bresenham(std::pair<int,int> start_coord, std::pair<int,int> end_coord){
        int x0 = start_coord.first, y0 = start_coord.second;
        int x1 = end_coord.first, y1 = end_coord.second;
        std::list<std::pair<int,int>> result;
        Orientation orientation(x1-x0,y1-y0);

        x1 -= x0;
        y1 -= y0;
        x1 = x1 < 0? -x1 : x1;
        y1 = y1 < 0? -y1 : y1;
        x0 = 0; y0 = 0;

        int dx =  x1 - x0, sx = x0 < x1 ? 1 : -1;
        int dy = -(y1 - y0), sy = y0 < y1 ? 1 : -1;

        // int dx =  abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
        // int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
        int err = dx + dy, e2;
        while (!(x0==x1 && y0==y1)){
            // std::cout << "(" << x0 << "," << y0 << ")" << std::endl;
            result.push_back(std::make_pair(x0,y0));
            e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
            if (e2 <= dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
        }

        orientation.transform(result);
        for (auto& i : result){
            i.first += start_coord.first;
            i.second += start_coord.second;
        }

        auto leading_it = result.begin(), following_it = result.begin();
        leading_it++;
        for (;leading_it!=result.end();++leading_it,++following_it){
            if ((leading_it->first-following_it->first)!=0 && (leading_it->second - following_it->first) != 0){
                result.insert(leading_it, std::make_pair(
                    leading_it->first - (leading_it->first - following_it->first),
                    leading_it->second
                    )
                );
                following_it++;
            }
        }

        return result;
        
    }

    static uint squared_distance_pixel(pix* a, pix* b){
        auto point_0 = get_xy(a->idx);
        auto point_1 = get_xy(b->idx);
        return  (uint) (point_1.first - point_0.first) * (point_1.first - point_0.first) +
                (point_1.second - point_0.second) * (point_1.second - point_0.second);
    }

    void set_condition(Condition _cond){
        cond = creator->create(_cond);
    }

    private: 
        static int nrows;
        static int ncols;
        static size_t numel;        
        static size_t num_valid_elements;    
        static pix* all_data;
        
        static std::unique_ptr<Creator> creator;  
        static std::unique_ptr<Up_condition> cond; 
            
};

int Data::nrows=0;
int Data::ncols=0;
size_t Data::numel=0;
size_t Data::num_valid_elements;
pix* Data::all_data = nullptr;
std::unique_ptr<Creator> Data::creator = std::make_unique<Creator>();
std::unique_ptr<Up_condition> Data::cond = Data::creator->create(Condition::Upstream);

// std::ostream &operator<<(std::ostream &os,const Data& m){
//         os << "elevation: " << m.elevation << " / flowacc: "  <<  m.flowacc << " / flowdir: " << m.flowdir << " / idx: " << m.idx <<  std::endl;
//         return os;
// }

   
struct Curtain{
    std::vector<pix*> clockwise;
    std::vector<pix*> counterclockwise;
    pix * center_pix;
    int squared_distance_clockwise, squared_distance_counterclockwise;
    void set_distances(){
        auto center = Data::get_xy(center_pix->idx);
        auto clockwise_coord = Data::get_xy((*(clockwise.rbegin()))->idx);
        auto counterclockwise_coord = Data::get_xy((*(counterclockwise.rbegin()))->idx);
        squared_distance_clockwise = (clockwise_coord.first - center.first) * (clockwise_coord.first - center.first) + 
                                     (clockwise_coord.second - center.second) * (clockwise_coord.second - center.second)  ;
        squared_distance_counterclockwise = (counterclockwise_coord.first - center.first) * (counterclockwise_coord.first - center.first) +
                                            (counterclockwise_coord.second - center.second) * (counterclockwise_coord.second - center.second);
    }
};


void bresenham(int x1, int y1, int x2, int y2)
{
    int m_new = 2 * (y2 - y1);
    int slope_error_new = m_new - (x2 - x1);
    for (int x = x1, y = y1; x <= x2; x++) {
        if ((x - x1) * (x - x1) + (y - y1)* (y - y1) <= 100 ) std::cout << "(" << x << "," << y << ")\n";
        else break;

        // Add slope to increment angle formed
        slope_error_new += m_new;

        // Slope error reached limit, time to
        // increment y and update slope error.
        if (slope_error_new >= 0) {
            y++;
            slope_error_new -= 2 * (x2 - x1);
        }
    }
}

int main(){

    // std::string tif_path="X:\\2023\\02\\20230265\\Calculations\\Hydrology\\Rasters\\Istad\\dem.tif";
    std::string tif_path="./data/dem.tif";
    std::cout << "Hell yeah!" << " " << tif_path <<" " << std::endl;

    //Registering all known GDAL drivers
    GDALAllRegister();
    CPLPushErrorHandler(CPLQuietErrorHandler);

    //Loading elevation, flow accumulation and flow direction   

    std::map<std::string,meta> all_meta;
    all_meta["elevation"] = meta("./data/dem.tif");
    all_meta["flow_dir"] = meta("./data/flowdir.tif");
    all_meta["flow_acc"] = meta("./data/flowacc.tif");

    // std::vector<Data> pixels;

    std::cout << "Getting raster info..." << std::endl;
    //Setting up storage:  this assumes that all the rasters are stored as float, have the same size, and have only one band
    GDALDataset *tif = (GDALDataset*) GDALOpen(all_meta["elevation"].file.c_str(), GA_ReadOnly );
    if(tif==NULL){
            std::cout << "Could not open file " << all_meta["elevation"].file << std::endl;
            return 1;
    }

    int nrows{tif->GetRasterYSize()},ncols{tif->GetRasterXSize()};
    size_t num_el =nrows * ncols;
    auto band = tif->GetRasterBand(1);
    auto nodata = (float) band->GetNoDataValue();

    float* const first_element = new float[num_el];

    auto cpl_err = band->RasterIO(GF_Read,0,0,ncols,nrows,(void*) first_element,ncols,nrows,tif-> GetRasterBand(1)->GetRasterDataType(),0,0);
    size_t num_data{0};
    float *pt = (float*) first_element;        
    for (size_t i = 0; i < num_el; ++i,pt++){
        if (*pt != nodata) num_data++;
    }

    std::cout << "Done!" << std::endl << "Loading data from rasters..." << std::endl;

    // Data* const pixels = new Data[num_data]; 
    Data data;
    Data::set_band_size(nrows,ncols);
    Data::set_number_valid_elements(num_data);
    // Data::set_pixels(pixels);

    
    
    pt = first_element;
    uint8_t * const save_data = new uint8_t[num_el];
    uint8_t* write_data = save_data;
    for (auto i=0; i<num_el;++i,++pt,++write_data){
        if (*pt ==  nodata){
            save_data[i]=-1;
        }
        else{
            save_data[i] = 0;
        }
    }

    //Trying to save to a tiff of bool
    auto driver_tif = GetGDALDriverManager()->GetDriverByName("GTiff");
    auto save_dataset =  driver_tif->Create("dummy.tif",ncols,nrows,1,GDT_Byte,NULL);
    double transform[6];
    tif->GetGeoTransform(transform);
    save_dataset->SetGeoTransform(transform);
    save_dataset->SetProjection( tif->GetProjectionRef() );
    save_dataset->GetRasterBand(1)->SetNoDataValue(-1);
    GDALClose(tif);
    
    

    pix* const all_data = new pix[num_data];
    Data::set_all_data(all_data);

    for (auto &[key,value] : all_meta){
        GDALDataset *tif;
        tif = (GDALDataset *) GDALOpen(value.file.c_str(), GA_ReadOnly);
        if(tif==NULL){
            std::cout << "Could not open file " << value.file << std::endl;
            return 1;
        }
        else std::cout << "Reading file " << value.file << std::endl;

        value.cols = tif -> GetRasterXSize();
        value.rows = tif -> GetRasterYSize();
        value.gdal_type = tif-> GetRasterBand(1)->GetRasterDataType();
        value.byte_size =  GDALGetDataTypeSizeBytes(value.gdal_type);

        band = tif->GetRasterBand(1);
        cpl_err = band->RasterIO(GF_Read,0,0,value.cols,value.rows,(void*) first_element,value.cols,value.rows,value.gdal_type,0,0);
        

        size_t cnt = 0;
        pt = first_element;

        if (key == "elevation"){
            for (size_t i = 0; i < num_el; ++i,pt++){
                if (*pt != nodata){
                    all_data[cnt].elevation = *pt;
                    all_data[cnt].idx = i;
                    cnt++;
                }
            }
        }
        else if (key == "flow_dir"){
            for (size_t i = 0; i < num_el; ++i,pt++){
                if (*pt != nodata){
                    all_data[cnt].dir = Direction{(uint)*pt};
                    cnt++;
                }
            }
        }
        else if (key == "flow_acc"){
            for (size_t i = 0; i < num_el; ++i,pt++){
                if (*pt != nodata){
                    all_data[cnt].flowacc = (size_t) *pt;
                    cnt++;   
                }
            }
        }
        GDALClose(tif);
    }

    delete[] first_element;
    std::cout << "Done!" << std::endl;




    auto start = std::chrono::high_resolution_clock::now();
    data.set_downstream_idx_for_all_data();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time setting index of downstream pixel in seconds: " << elapsed.count() << std::endl;


    start = std::chrono::high_resolution_clock::now();
    Data::set_all_upstream_pixels();
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Elapsed time setting all upstream pixels in seconds: " << elapsed.count() << std::endl;

    
    start = std::chrono::high_resolution_clock::now();
    Data::set_upstream_with_largest_flowacc();
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Elapsed time finding upstream pixel with largest flow accumulation in seconds: " << elapsed.count() << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    // size_t idx = Data::get_idx(940,9310);
    size_t idx = Data::get_idx(2286,7842);
    pix* outlet = Data::binary_search_pix(idx);
    data.set_condition(Condition::Dam);
    data.get_basin(outlet);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Elapsed time setting upstream pixels in seconds: " << elapsed.count() << std::endl;
    data.paint_pixels(save_data);

    start = std::chrono::high_resolution_clock::now();
    data.set_condition(Condition::River);
    data.get_basin(outlet);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Elapsed time setting upstream pixels in seconds: " << elapsed.count() << std::endl;
    data.paint_pixels(save_data);

    
    // auto start_coord = Data::get_xy((*(data.basin.begin()))->idx);
    // auto end_coord =  Data::get_xy((*(data.basin.rbegin()))->idx);

    // std::cout << start_coord.first << " " << start_coord.second << std::endl << end_coord.first << " " << end_coord.second << std::endl;
    // std::cout << end_coord.first - start_coord.first << " " << end_coord.second - start_coord.second << std::endl;
    // bresenham(0,0,20,3);

    // bresenham(start_coord.first,start_coord.second,end_coord.first,end_coord.second);




    auto curtain = data.get_curtain(*(data.basin.begin()),*(data.basin.rbegin()));

    for (auto& p: curtain){
        save_data[p->idx]=1;
    }

    auto save_err = save_dataset->GetRasterBand(1)->RasterIO(GF_Write,0,0,ncols,nrows,
                                                            (void*) save_data,ncols,nrows
                                                            ,GDT_Byte,0,0);

    
    GDALClose(save_dataset);
    delete[] save_data;

    // std::cout << "I: " << std::endl;
    // bresenham(0,0,10,10);
    // std::cout << "II: " << std::endl;
    // bresenham(0,0,10,-10);
    // std::cout << "III: " << std::endl;
    // bresenham(0,0,-10,-10);
    // std::cout << "IV: " << std::endl;
    // bresenham(0,0,-10,10);

    // size_t cnt = 0;
    // for (auto& i : pixels){
    //     if (i.flowacc > 1000) ++cnt;
    // } 
    // std::cout << "There are " << cnt << " pixels with the adequate conditions" << std::endl;

    // std::cout << "Are the pixels sorted: " << std::boolalpha << std::is_sorted(pixels.begin(),pixels.end()) << std::endl;
    

    // delete [] pixels;
    delete [] all_data;
    // delete[] first_element;
    return 0;
} 