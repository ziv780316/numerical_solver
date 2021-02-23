###QA_TAG matrix
###QA_GOLDEN_FILE matrix

# paramters
cp -f ${DEVELOP_ROOT_DIR}/numerical_solver/matrix_solver/regression/Makefile .
make
./test > matrix

