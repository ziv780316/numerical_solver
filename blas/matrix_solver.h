#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

// y := alpha*A*x + beta*y, or y := alpha*A*x + beta*y if transpose
int dense_matrix_vector_multiply ( int m, int n, double alpha, double *A, double *x, double beta, double *y, bool transpose );

// C := alpha*A*B + beta*C,
int dense_matrix_matrix_multiply ( int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, bool a_transpose, bool b_transpose );

// solve A*x = b, or (A**T)*x = b where A is triangular matrix, x is RHS and result will overwrite in x after solve
int dense_triangular_solve ( int n, double *A, double *x, bool is_lower_triangular, bool transpose, bool is_unit_triangular );

// print matrix
int dense_print_matrix ( int m, int n, double *A );

#endif

