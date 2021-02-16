#!/bin/bash
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/suitesparse/lib
make -j 4 LDFLAGS+='-L${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/suitesparse/lib'   LAPACK='-L${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib -llapack -lgfortran' BLAS='-L${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib -lblas -lgfortran'
