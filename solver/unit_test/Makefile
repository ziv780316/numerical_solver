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
BUILD_TYPE = -O3
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -std=gnu99 -fPIC 
LIBS       = -L../lib -lmatrix_solver -lblas -llapack -lgfortran -lm
INCLUDE    = -I.. -I.
MOVE       = mv -f

#=================================================
# Build target
#=================================================
LIBRARY_DIR       := library
UNIT_TEST_OBJ     := dense_matrix_vector_multiply.o dense_matrix_matrix_multiply.o dense_triangular_solve.o dense_swap_vector.o dense_lu_factor.o dense_factor_and_solve.o dense_matrix_scale.o dense_matrix_rank_1_update.o dense_matrix_inverse.o

#=================================================
# Compile implicit rules
#=================================================
%.o:%.c
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) 
	$(LINKER) $@ -o $(@:.o=) $(LIBS)
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command


#=================================================
# Target rules
#=================================================
default: unit_test

all: unit_test 

debug: BUILD_TYPE:=-g3 
debug: default

create_dir:
	@mkdir -p $(LIBRARY_DIR) 

unit_test: create_dir $(UNIT_TEST_OBJ)

clean:
	rm -rf $(LIBRARY_DIR)/*


