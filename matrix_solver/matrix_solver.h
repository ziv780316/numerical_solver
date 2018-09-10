#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

typedef enum
{
	REAL_NUMBER,
	COMPLEX_NUMBER
} number_type;

typedef enum
{
	TRANS_NONE,
	TRANS_NORMAL,
	TRANS_CONJUGATE
} transpose_type;

typedef struct
{
	double real;
	double imag;
} complex_t;

// x := alpha*x
int dense_vector_scale ( int n, double *x, double *alpha, number_type );

// aij = aii * alpha
int dense_matrix_scale ( int m, int n, double *A, double *alpha, number_type );

// obtain D of A=L+D+U 
int dense_matrix_get_diagonal ( int n, double *A, double *D, number_type );

// main diagonal addition, aii = aii + alpha
int dense_diagonal_addition ( int n, double *A, double *alpha, number_type );

// val = x . y
int dense_vector_inner_product ( int n, double *x, double *y, double *val, bool conjugate, number_type type );

// A = x . yT
int dense_vector_outer_product ( int m, int n, double *x, double *y, double *A, bool conjugate, number_type );

// A = A + alpha * (x . yT)
int dense_maxtrix_rank_1_update ( int m, int n, double *A, double *alpha, double *x, double *y, bool conjugate, number_type );

// vector norm |x|p = (sum |xi|^p)^(1.0/p)
int dense_vector_norm ( int p_norm, int n, double *x, double *val, number_type );

// y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y or y := alpha*A**H*x + beta*y
int dense_matrix_vector_multiply ( int m, int n, double *alpha, double *A, double *x, double *beta, double *y, transpose_type , number_type );

// C := alpha*A*B + beta*C,
int dense_matrix_matrix_multiply ( int m, int n, int k, double *alpha, double *A, double *B, double *beta, double *C, bool a_transpose, bool b_transpose, number_type );

// solve A*x = b, or (A**T)*x = b where A is triangular matrix, x is RHS and result will overwrite in x after solve
int dense_triangular_solve ( int n, double *A, double *x, bool is_lower_triangular, bool transpose, bool is_unit_triangular, number_type );

// swap vector x and y
int dense_swap_vector ( int n, double *x, double *y, number_type );

// LU factorization A = P * L * U
int dense_lu_factor ( int n, double *A, int *p, number_type );

// solve A*x = b 
int dense_solve ( int n, double *A, double *x, int *p, bool transpose, number_type );

// solve A*x = b with LU factorization
int dense_factor_and_solve ( int n, double *A, double *x, bool transpose, number_type );

// A := A^-1, computes inv(A) by solving the system inv(A)*L = inv(U) for inv(A).
int dense_matrix_inverse ( int n, double *A, int *p, number_type );

// print matrix
int dense_print_vector ( int n, double *x, number_type );
int dense_print_vector_i ( int n, int *x, number_type );
int dense_print_matrix ( int m, int n, double *A, number_type );
int dense_print_matrix_LU ( int n, double *A, number_type );

#endif

