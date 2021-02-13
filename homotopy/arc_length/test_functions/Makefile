#=================================================
# Makefile settings
#=================================================
.PHONY: clean # specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=../library # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
LINKER     = $(CC)
DEFINES    = 
BUILD_TYPE = -O3
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -Wextra -std=gnu99 -fPIC 
MAT_SOLVER_PATH = ${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib
MAT_SOLVER_INC = ${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver
LIBS       = -lm -L$(MAT_SOLVER_PATH) -lmatrix_solver ${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib/libblas.so ${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib/liblapack.so -lgfortran -ldl 
INCLUDE    = -I. -I$(MAT_SOLVER_INC)
MOVE       = mv -f

#=================================================
# Build target
#=================================================
TEST_FUNCTION_OBJ := f1.so f2.so f3.so f4.so f5.so f6.so f7.so f8.so f9.so

#=================================================
# Compile implicit rules
#=================================================
%.so:%.c
	$(CC) -shared $< -o $@ $(CFLAGS) $(INCLUDE) ${LIBS}


#=================================================
# Target rules
#=================================================
default: $(TEST_FUNCTION_OBJ)

debug: BUILD_TYPE:=-g3 
debug: default

clean:
	rm -rf $(TEST_FUNCTION_OBJ)


