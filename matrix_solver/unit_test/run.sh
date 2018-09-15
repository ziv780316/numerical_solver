#!/bin/bash

functions="\
dense_factor_and_solve \
dense_lu_factor \
dense_matrix_inverse \
dense_matrix_matrix_multiply \
dense_matrix_rank_1_update \
dense_matrix_scale \
dense_matrix_vector_multiply \
dense_swap_vector \
dense_triangular_solve \
complex_dense_vector_scale \
complex_dense_matrix_scale \
complex_dense_vector_inner_product \
complex_dense_matrix_rank_1_update \
complex_dense_vector_norm \
complex_dense_matrix_vector_multiply \
complex_dense_matrix_matrix_multiply \
complex_dense_triangular_solve \
complex_dense_swap_vector \
complex_dense_lu_factor \
"
 
# execute each unit function
for f in $functions; do
	./bin/$f > ${f}.log
done

# compare with golden
for f in $functions; do
	diff -q ${f}.log golden/${f}.log
done

# memory test
#dynamic_analyis="valgrind --leak-check=full --show-leak-kinds=all -v"
#for f in $functions; do
#	$dynamic_analyis ./bin/$f >& ${f}.valgrind
#done
