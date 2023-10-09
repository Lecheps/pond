#include <iostream>
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

struct pix{
    float elevation;
    Direction dir; 
    size_t flowacc;
    size_t idx;
    pix* downstream;
    std::vector<pix*>upstream;
    pix* max_upstream;
};

enum Condition {Upstream,River,Dam};     

class Up_condition{
    public:
        virtual bool check(const pix* up, const pix* down) = 0;
        // virtual void update_condition(float) = 0;
        virtual ~Up_condition() = default; 
        // static float outlet_elevation;
        // static float highest_surrounding_elevation;
        static bool follow_steepest_path;
        static pix* outlet;
};

// float Up_condition::outlet_elevation = 0.0;
// float Up_condition::highest_surrounding_elevation;
bool Up_condition::follow_steepest_path = false;
pix* Up_condition::outlet = nullptr;

class Up_condition_upstream: public Up_condition {
    public:
        bool check(const pix* up, const pix* down){
            return (up->downstream == down);
        }
        // void update_condition(float val) {};
        
        
};

class Up_condition_dam: public Up_condition {
    public:
        bool check(const pix* up, const pix* down){
            return (up->downstream == down && (up->elevation - outlet->elevation) <= 2.0);
        }
        // void update_condition(float val) {outlet_elevation = val;};
        
}; 

class Up_condition_river: public Up_condition {
    public:
        bool check(const pix* up, const pix* down){
            return (true);
        }
        // void update_condition(float val) {outlet_elevation = val;};
        Up_condition_river() {follow_steepest_path=true;};
};

class Creator {
    public:
        virtual std::unique_ptr<Up_condition> create(Condition cond){
            if(Condition::Upstream == cond) return std::make_unique<Up_condition_upstream>();
            if(Condition::Dam == cond) return std::make_unique<Up_condition_dam>();
            if(Condition::River == cond) return std::make_unique<Up_condition_river>();
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
            // cond->update_condition(outlet->elevation);
            cond->outlet = outlet;
            get_all_upstream(outlet);
            std::cout << "Number of upstream pixels for this pixel: " << basin.size() << std::endl;
        }

        void get_all_upstream(pix* data){
            // float max_surrounding_elev  = std::numeric_limits<float>::min();
            // for (auto &i: data->upstream){
            //     if (i->elevation > max_surrounding_elev) max_surrounding_elev = i->elevation; 
            // }
            // cond->update_condition(max_surrounding_elev);
            if (!cond->follow_steepest_path){
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
                std::cout << data->idx << " " << data->flowacc << std::endl;
                if (data->max_upstream != nullptr) get_all_upstream(data->max_upstream);
            }
                      
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
    


    private: 
        static int nrows;
        static int ncols;
        static size_t numel;        
        static size_t num_valid_elements;    
        static pix* all_data;
        std::vector<pix*> basin;
        std::unique_ptr<Creator> creator = std::make_unique<Creator>();
        std::unique_ptr<Up_condition> cond = creator->create(Condition::River);
            
};

int Data::nrows=0;
int Data::ncols=0;
size_t Data::numel=0;
size_t Data::num_valid_elements;
pix* Data::all_data = nullptr;

// std::ostream &operator<<(std::ostream &os,const Data& m){
//         os << "elevation: " << m.elevation << " / flowacc: "  <<  m.flowacc << " / flowdir: " << m.flowdir << " / idx: " << m.idx <<  std::endl;
//         return os;
// }





 

   

    void bresenham(int x1, int y1, int x2, int y2)
    {
        int m_new = 2 * (y2 - y1);
        int slope_error_new = m_new - (x2 - x1);
        for (int x = x1, y = y1; x <= x2; x++) {
            std::cout << "(" << x << "," << y << ")\n";
    
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

    GDALClose(tif);
    std::cout << "Done!" << std::endl << "Loading data from rasters..." << std::endl;

    // Data* const pixels = new Data[num_data]; 
    Data data;
    Data::set_band_size(nrows,ncols);
    Data::set_number_valid_elements(num_data);
    // Data::set_pixels(pixels);
    

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
                    // pixels[cnt].elevation = *pt;
                    // pixels[cnt].idx=i;
                    all_data[cnt].idx = i;
                    cnt++;
                }
            }
        }
        else if (key == "flow_dir"){
            for (size_t i = 0; i < num_el; ++i,pt++){
                if (*pt != nodata){
                    all_data[cnt].dir = Direction{(uint)*pt};
                    // pixels[cnt].flowdir = (uint) *pt;
                    // pixels[cnt].idx=i;
                    cnt++;
                }
            }
        }
        else if (key == "flow_acc"){
            for (size_t i = 0; i < num_el; ++i,pt++){
                if (*pt != nodata){
                    all_data[cnt].flowacc = (size_t) *pt;
                    // pixels[cnt].flowacc = (size_t) *pt;
                    // pixels[cnt].idx=i;
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
    for (auto i = 0; i < num_data; ++i){
        pix* pt;
        // float max_elev = std::numeric_limits<float>::min();
        size_t max_accumulation = 0;
        if (! all_data[i].upstream.empty()) {
            for (std::vector<pix*>::iterator j = all_data[i].upstream.begin(); j !=all_data[i].upstream.end(); ++j) {
                // if ((*j)->elevation > max_elev) pt = *j;
                if ((*j)->flowacc > max_accumulation) {
                    pt = *j;
                    max_accumulation = (*j)->flowacc;
                }
            }
        } else pt = nullptr;
        
        all_data[i].max_upstream = pt;
        
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Elapsed time finding steepest path in seconds: " << elapsed.count() << std::endl;
    
    
    
    // for (size_t i = 0; i < 100; ++i){
    //     std::cout << all_data[i].idx  << std::endl << "\t";
    //     for (auto &j : all_data[i].upstream){
    //         std::cout << j->idx << " ";
    //     }
    //     std::cout << std::endl;
    // }

    start = std::chrono::high_resolution_clock::now();
    // size_t idx = Data::get_idx(38,9662);
    size_t idx = Data::get_idx(940,9310);
    pix* outlet = Data::binary_search_pix(idx);
    data.get_basin(outlet);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Elapsed time setting upstream pixels in seconds: " << elapsed.count() << std::endl;


    bresenham(0,0,20,3);

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