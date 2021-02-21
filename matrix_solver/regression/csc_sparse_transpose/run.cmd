###QA_TAG matrix
###QA_GOLDEN_FILE matrix.out

# paramters
cp -f ../Makefile .
make
./test > matrix

