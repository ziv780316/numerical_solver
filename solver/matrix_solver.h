#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

// x := alpha*x
int dense_vector_scale ( int n, double *x, double alpha );

// aij = aii * alpha
int dense_matrix_scale ( int m, int n, double *A, double alpha );

// main diagonal addition, aii = aii + alpha
int dense_diagonal_addition ( int n, double *A, double alpha );

// val = x . y
int dense_vector_inner_product ( int n, double *x, double *y, double *val );

// A = x . yT
int dense_vector_outer_product ( int n, double *x, double *y, double *A );

// A = A + alpha * (x . yT)
int dense_maxtrix_rank_1_update ( int n, double *A, double alpha, double *x, double *y );

// vector norm |x|p = (sum |xi|^p)^(1.0/p)
int dense_vector_norm ( int p_norm, int n, double *x, double *val );

// y := alpha*A*x + beta*y, or y := alpha*AT*x + beta*y if transpose
int dense_matrix_vector_multiply ( int m, int n, double alpha, double *A, double *x, double beta, double *y, bool transpose );

// C := alpha*A*B + beta*C,
int dense_matrix_matrix_multiply ( int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, bool a_transpose, bool b_transpose );

// solve A*x = b, or (A**T)*x = b where A is triangular matrix, x is RHS and result will overwrite in x after solve
int dense_triangular_solve ( int n, double *A, double *x, bool is_lower_triangular, bool transpose, bool is_unit_triangular );

// swap vector x and y
int dense_swap_vector ( int n, double *x, double *y );

// LU factorization A = P * L * U
int dense_lu_factor ( int n, double *A, int *p );

// solve A*x = b 
int dense_solve ( int n, double *A, double *x, int *p, bool transpose );

// solve A*x = b with LU factorization
int dense_factor_and_solve ( int n, double *A, double *x, bool transpose );

// A := A^-1, computes inv(A) by solving the system inv(A)*L = inv(U) for inv(A).
int dense_matrix_inverse ( int n, double *A, int *p );

// print matrix
int dense_print_vector ( int n, double *x );
int dense_print_vector_i ( int n, int *x );
int dense_print_matrix ( int m, int n, double *A );
int dense_print_matrix_LU ( int n, double *A );

#endif

