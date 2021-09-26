#=================================================
# Makefile settings
#=================================================
.PHONY: create_dir # specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=library # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
LINKER     = $(CC)
DEFINES    = 
BUILD_TYPE = -g3 -g
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -Wextra -Wenum-compare -std=gnu99 -fPIC -Wno-unused-parameter -Wno-unused-variable
FFTW_LIBS  = -Wl,-rpath '/usr/lib/x86_64-linux-gnu' -lfftw3
FFTW_INC_PATH = /usr/include
PY_LIBS       = /home/ziv/work1/compiler/c_python/library/c_python.so -Wl,-rpath '/lib/x86_64-linux-gnu' -lpython3.8 
PY_INCLUDE    = -I/usr/include/python3.8 -I/home/ziv/work1/compiler/c_python
LIBS       = -lm $(FFTW_LIBS) $(PY_LIBS)
INCLUDE    = -I$(FFTW_INC_PATH) $(PY_INCLUDE)
MOVE       = mv -f

#=================================================
# Build target
#=================================================
LIBRARY_DIR       := library
TEST_FUNCTION_DIR := test_functions
BIN_DIR       	  := bin
OBJS              := main.o opts.o conv.o 
TARGET            := convolution_solver gen_response

#=================================================
# Compile implicit rules
#=================================================
%.o:%.c
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command

%.e:%.c
	$(CC) -E $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR)

%.o:%.cxx
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR)

%.e:%.cxx
	$(CXX) -E $< -o $@ $(CXXFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR)


#=================================================
# Target rules
#=================================================
default: $(TARGET)

convolution_solver: create_dir $(OBJS) 
	cd $(LIBRARY_DIR);\
	$(LINKER) -o $@ $(OBJS) $(LIBS);\
	$(MOVE) $@ ../$(BIN_DIR)

gen_response: gen_response.o
	$(LINKER) $(LIBRARY_DIR)/gen_response.o -o $@;\
	$(MOVE) $@ $(BIN_DIR)

debug: BUILD_TYPE:=-g3 
debug: default

create_dir:
	@mkdir -p $(LIBRARY_DIR)
	@mkdir -p $(BIN_DIR)

unit_test: create_dir $(UNIT_TEST_OBJ)

clean:
	rm -rf $(LIBRARY_DIR)/*


