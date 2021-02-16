#!/bin/bash
make -j 4 LAPACK=${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib/liblapack.so BLAS=${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/lib/libblas.so
