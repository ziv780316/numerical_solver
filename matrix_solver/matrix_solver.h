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
	TRANS_CONJUGATE // Hermitian transpose 
} transpose_type;

typedef enum
{
	TRIG_LOWER,
	TRIG_UPPER,
	TRIG_LOWER_UNIT, // main diagonal is 1
	TRIG_UPPER_UNIT
} triangular_type;

typedef enum
{
	FACTOR_LU_CROUT,
	FACTOR_LU_SIVAN_RECURSIVE,
	FACTOR_LU_LEFT_LOOKING,
	FACTOR_LU_RIGHT_LOOKING,
	FACTOR_LU_CHOLESKY, // A = L*(L**H), need real symmetric (hermitian) positive definite matrix
	FACTOR_LU_BUNCH_KAUFMAN // A = L*D*(L**T), D is block diagonal matrix, need real symmetric (hermitian) matrix 
} factorization_type;

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

// C := alpha*OP(A)*OP(B) + beta*C, OP(A) = A or A**T or A**H, OP(A)*OP(B) is m*k*n
int dense_matrix_matrix_multiply ( int ma, int na, int mb, int nb, double *alpha, double *A, double *B, double *beta, double *C, transpose_type a_transpose, transpose_type b_transpose, number_type type );

// solve A*x = b, or (A**T)*x = b where A is triangular matrix, x is RHS and result will overwrite in x after solve
int dense_triangular_solve ( int n, double *A, double *x, triangular_type, transpose_type, number_type );

// swap vector x and y
int dense_swap_vector ( int n, double *x, double *y, number_type );

// LU factorization A = P * L * U = Pinv^-1 * L * U = Pinv' * L * U
int dense_lu_factor ( int n, double *A, int *pinv, factorization_type, number_type );

// solve A*x = b 
int dense_solve ( int n, int nrhs, double *A, double *x, int *pinv, factorization_type, transpose_type transpose, number_type );

// eval |A|
double dense_eval_det ( int n, double *A, int *pinv, factorization_type, number_type );

// solve A*x = b with LU factorization
int dense_factor_and_solve ( int n, int nrhs, double *A, double *x, int *pinv, factorization_type, transpose_type, number_type );

// A := A^-1, computes inv(A) by solving the system inv(A)*L = inv(U) for inv(A).
int dense_matrix_inverse ( int n, double *A, int *pinv, factorization_type, number_type );

// |A|_1 = col sum, |A|_max = row sum ...
int dense_matrix_norm ( int p_norm, int m, int n, double *A, double *val, number_type );

// print matrix
int dense_print_vector ( int n, double *x, number_type );
int dense_print_vector_i ( int n, int *x );
int dense_print_matrix ( int m, int n, double *A, number_type );
int dense_print_matrix_perm ( int n, int *p );
int dense_print_matrix_LU ( int n, double *A, number_type );
int dense_print_matrix_trig ( int n, double *A, triangular_type, number_type );

extern int g_matrix_print_format;
#define MATRIX_PRINT_FORMAT_PLAIN 0
#define MATRIX_PRINT_FORMAT_MATLAB 1

#endif

