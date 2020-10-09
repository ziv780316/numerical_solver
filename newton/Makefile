#=================================================
# Makefile settings
#=================================================
.PHONY: create_dir # specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

#VPATH=library # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
LINKER     = $(CC)
DEFINES    = 
BUILD_TYPE = -g3
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -Wextra -Wenum-compare -std=gnu99 -fPIC -Wno-unused-parameter -Wno-unused-variable
MAT_SOLVER_PATH = /home/lintu/baba/numerical/matrix_solver/lib
LIBS       = -L$(MAT_SOLVER_PATH) -lmatrix_solver -llapack -lblas -lgfortran -lm -ldl
INCLUDE    = -I../matrix_solver -I. 
MOVE       = mv -f

#=================================================
# Build target
#=================================================
LIBRARY_DIR       := library
TEST_FUNCTION_DIR := test_functions
BIN_DIR       	  := bin
OBJS              := newton.o main.o opts.o
TARGET            := newton_solver libnewton.a

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

libnewton.a: $(LIBRARY_DIR)/newton.o
	ar -rv $@ $^
	mv $@ lib

newton_solver: create_dir $(OBJS) make_test_functions
	cd $(LIBRARY_DIR);\
	$(LINKER) -o $@ $(OBJS) $(LIBS);\
	$(MOVE) $@ ../$(BIN_DIR)

make_test_functions: 
	$(MAKE) -C $(TEST_FUNCTION_DIR) OPTS=$(BUILD_TYPE) 

debug: BUILD_TYPE:=-g3 
debug: default

create_dir:
	@mkdir -p lib
	@mkdir -p $(LIBRARY_DIR)
	@mkdir -p $(BIN_DIR)

unit_test: create_dir $(UNIT_TEST_OBJ)

clean:
	rm -rf $(LIBRARY_DIR)/*
	rm -rf $(BIN_DIR)/*


