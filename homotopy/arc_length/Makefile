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
BUILD_TYPE = -O2 -g
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -Wextra -Wenum-compare -std=gnu99 -fPIC -Wno-unused-parameter -Wno-unused-variable
MAT_SOLVER_PATH = /home/lintu/baba/numerical/matrix_solver/lib
NEWTON_SOLVER_PATH = /home/lintu/baba/numerical/newton/lib
NEWTON_SOLVER_INC = /home/lintu/baba/numerical/newton
MAT_SOLVER_INC = /home/lintu/baba/numerical/matrix_solver
LIBS       = -L$(NEWTON_SOLVER_PATH) -lnewton -L$(MAT_SOLVER_PATH) -lmatrix_solver -llapack -lblas -lgfortran -lm -ldl 
INCLUDE    = -I$(MAT_SOLVER_INC) -I$(NEWTON_SOLVER_INC) -I. 
MOVE       = mv -f

#=================================================
# Build target
#=================================================
LIBRARY_DIR       := library
TEST_FUNCTION_DIR := test_functions
BIN_DIR       	  := bin
OBJS              := main.o arc_length.o opts.o
TARGET            := homotopy_solver 

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

homotopy_solver: create_dir $(OBJS) make_test_functions
	cd $(LIBRARY_DIR);\
	$(LINKER) -o $@ $(OBJS) $(LIBS);\
	$(MOVE) $@ ../$(BIN_DIR)

make_test_functions: 
	$(MAKE) -C $(TEST_FUNCTION_DIR) OPTS="$(BUILD_TYPE)"

debug: BUILD_TYPE:=-g3 
debug: default

create_dir:
	@mkdir -p lib
	@mkdir -p $(LIBRARY_DIR)
	@mkdir -p $(BIN_DIR)

unit_test: create_dir $(UNIT_TEST_OBJ)

clean:
	rm -rf $(LIBRARY_DIR)/*


