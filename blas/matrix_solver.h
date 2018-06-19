#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

// y := alpha*A*x + beta*y, or y := alpha*(A**T)*x + beta*y if transpose
int dense_matrix_vector_multiply ( int m, int n, double alpha, double *A, double *x, double beta, double *y, bool transpose );

#endif

