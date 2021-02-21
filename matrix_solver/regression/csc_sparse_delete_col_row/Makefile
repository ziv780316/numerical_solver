#=================================================
# Makefile settings
#=================================================
.DEFAULT_GOAL=default # specify default target if there is no target specified

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
LINKER     = $(CC)
DEFINES    = 
BUILD_TYPE = -g3
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -std=gnu99 -fPIC 
LIBS       = -L${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib -lmatrix_solver -llapack -lblas -lgfortran -lm -L${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/suitesparse/lib -lcxsparse 
INCLUDE    = -I${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver -I${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/suitesparse/include
MOVE       = mv -f

#=================================================
# Build target
#=================================================
UNIT_TEST_CASE = test.o ${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib/libmatrix_solver.a

#=================================================
# Compile implicit rules
#=================================================
%.o:%.c 
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) 
	$(LINKER) $@ -o $(@:.o=) $(LIBS) 


#=================================================
# Target rules
#=================================================
default: $(UNIT_TEST_CASE)

debug: BUILD_TYPE:=-g3 
debug: default

clean:
	rm -rf $(UNIT_TEST_CASE) $(UNIT_TEST_CASE:.o=) $(UNIT_TEST_CASE:.o=.c)


