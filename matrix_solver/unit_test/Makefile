#=================================================
# Makefile settings
#=================================================
.PHONY: tags create_dir # specify phony target, do not generate file as target name

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
MATRIX_SOLVER_LIB = ../lib/libmatrix_solver.a
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -std=gnu99 -fPIC 
LIBS       = -L../lib -lmatrix_solver -llapack -lblas -lgfortran -lm
INCLUDE    = -I.. -I.
MOVE       = mv -f

#=================================================
# Build target
#=================================================
BIN_DIR       	  := bin
LIBRARY_DIR       := library
UNIT_TEST_OBJ     := dense_matrix_vector_multiply.o\
		     dense_matrix_matrix_multiply.o\
		     dense_triangular_solve.o\
		     dense_swap_vector.o\
		     dense_lu_factor.o\
		     dense_factor_and_solve.o\
		     dense_matrix_scale.o\
		     dense_matrix_rank_1_update.o\
		     dense_matrix_inverse.o\
		     complex_dense_vector_scale.o\
		     complex_dense_matrix_scale.o\
		     complex_dense_vector_inner_product.o\
		     complex_dense_matrix_rank_1_update.o\
		     complex_dense_vector_norm.o\
		     complex_dense_matrix_vector_multiply.o\
		     complex_dense_matrix_matrix_multiply.o\
		     complex_dense_triangular_solve.o\
		     complex_dense_swap_vector.o\
		     complex_dense_lu_factor.o\

#=================================================
# Compile implicit rules
#=================================================
%.o:%.c $(MATRIX_SOLVER_LIB) 
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) 
	$(LINKER) $@ -o $(@:.o=) $(LIBS) 
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command
	@$(MOVE) $(@:.o=) $(BIN_DIR)


#=================================================
# Target rules
#=================================================
default: unit_test

all: unit_test 

debug: BUILD_TYPE:=-g3 
debug: default

create_dir:
	@mkdir -p $(BIN_DIR) 
	@mkdir -p $(LIBRARY_DIR) 

unit_test: create_dir $(UNIT_TEST_OBJ)

tags: 
	rm -f tags
	find ../$(BLAS_VERSION) -name "*.f" | xargs ctags --language-force=fortran -f /dev/stdout | sed -r 's/([a-z]\S+)(\s+.*(subroutine|function))/\L\1\E\2/I' >> tags
	find ../$(LAPACK_VERSION) -name "*.f" | xargs ctags --language-force=fortran -f /dev/stdout | sed -r 's/([a-z]\S+)(\s+.*(subroutine|function))/\L\1\E\2/I' >> tags
	find .. -name "*.c" -o -name "*.h" | xargs ctags --append 

clean:
	rm -rf $(BIN_DIR)/*
	rm -rf $(LIBRARY_DIR)/*


