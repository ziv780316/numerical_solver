#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

// y := alpha*A*x + beta*y, or y := alpha*A*x + beta*y if transpose
int dense_matrix_vector_multiply ( int m, int n, double alpha, double *A, double *x, double beta, double *y, bool transpose );

// C := alpha*A*B + beta*C,
int dense_matrix_matrix_multiply ( int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, bool a_transpose, bool b_transpose );

#endif

