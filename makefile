#IDIR=/usr/include/gdal
#CC=/usr/bin/g++
#EXE=ponds

#DBGDIR=debug
#DBGEXE=$(DBGDIR)/$(EXE)
#DBGOBJS=$(addprefix $(DBGDIR)/, $(OBJS))
#DBGCFLAGS=-g -O0 -DDEBUG

#RELDIR=release
#RELEXE=$(RELDIR)/$(EXE)
#RELOBJS=$(addprefix $(RELDIR)/, $(OBJS))
#RELCFLAGS=-g -O3 -DNDEBUG


#.PHONY: all clean debug prep release remake

#all: prep release

#debug: $(DBGEXE)

#$(DBGEXE): $(DBGOBJS)
#    $(CC) $(CFLAGS) $(DBGCFLAGS) -o $(DBGEXE) $^

#$(DBGDIR)/%.o: %.cpp
#    $(CC) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $<

#
# Release rules
#
#release: $(RELEXE)

#$(RELEXE): $(RELOBJS)
#    $(CC) $(CFLAGS) $(RELCFLAGS) -o $(RELEXE) $^

#$(RELDIR)/%.o: %.c
#    $(CC) -c $(CFLAGS) $(RELCFLAGS) -o $@ $<

#
# Other rules
#
#prep:
#    @mkdir -p $(DBGDIR) $(RELDIR)

#remake: clean all

#clean:
#    rm -f $(RELEXE) $(RELOBJS) $(DBGEXE) $(DBGOBJS)

SRC_DIR:=src
OBJ_DIR:=obj
BIN_DIR:=bin

EXE:=$(BIN_DIR)/pond
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o) main.o

CC=g++
CPPFLAGS:=-Iinclude -I/usr/include/gdal -MMD -P
CFLAGS:=-Wall -g
LDFLAGS:=-L/usr/lib
LDLIBS=-lgdal

.PHONY: all clean

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp main.cpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR) main.o main.d

-include $(OBJ:.o=.d)


