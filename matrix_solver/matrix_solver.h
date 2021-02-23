#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

#include  "cs.h"
#include  "klu.h"

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

// ----------------------------------------------------------------------
// Dense Matrix Operation (Column Major, 0-base)
// ----------------------------------------------------------------------

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
double dense_eval_factor_det ( int n, double *A, int *pinv, factorization_type, number_type );
double dense_eval_det ( int n, double *A, factorization_type, number_type );

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

// ----------------------------------------------------------------------
// Sparse Matrix Operation (CSC, 0-base)
// ----------------------------------------------------------------------
typedef long sparse_int; 
typedef double sparse_float; 
typedef struct
{
	sparse_int row;
	sparse_int col;
	sparse_float x;
} sparse_element_t;

typedef struct
{
	sparse_int nz;
	sparse_int m;
	sparse_int n;
	sparse_element_t *elements;
} sparse_triplet_t;

typedef struct
{
	sparse_int nz;
	sparse_int m;
	sparse_int n;
	sparse_int *Ap;
	sparse_int *Ai;
	sparse_float *Ax;
} sparse_csc_t;

typedef enum
{
	SPARSE_LU_DECOMPOSITION_KLU,
} sparse_lu_method;

// basic
void delete_sparse ( sparse_csc_t *A );
void alloc_sparse ( sparse_csc_t *A );
sparse_csc_t *copy_sparse ( sparse_csc_t *A );
void copy_csc_to_CXSparseCSC ( sparse_csc_t *A, cs_dl *B );
void copy_CXSparseCSC_to_csc ( cs_dl *A, sparse_csc_t *B );
sparse_csc_t *sparse_convert_triplet_to_CSC ( sparse_triplet_t *A );

// scale, A = A*alpha
void sparse_matrix_scale ( sparse_csc_t *A, sparse_float alpha );

// addition, C = alpha*A + beta*B
sparse_csc_t *sparse_matrix_addition ( sparse_csc_t *A, sparse_csc_t *B, sparse_float alpha, sparse_float beta );

// addition, C = A * B
sparse_csc_t *sparse_matrix_matrix_multiply ( sparse_csc_t *A, sparse_csc_t *B );

// Ax = b
void sparse_matrix_multiply_vector ( sparse_csc_t *A, sparse_float *x, sparse_float *b );

// A = inv(P)*R*L*U*Q
int sparse_matrix_lu_decomposition ( sparse_csc_t *A, sparse_lu_method method, sparse_csc_t *L, sparse_csc_t *U, sparse_csc_t *Pinv, sparse_csc_t *Q, sparse_csc_t *R );

// tranpose 
int sparse_matrix_transpose ( sparse_csc_t *A );

// convert CSC A to CSR A
int sparse_matrix_csc_to_csr ( sparse_csc_t *A );

// delete row or col 
int sparse_matrix_delete_col ( sparse_csc_t *A, sparse_int col );
int sparse_matrix_delete_row ( sparse_csc_t *A, sparse_int row );

// print matrix
sparse_float *sparse_to_full_matrix ( sparse_csc_t *A );
int sparse_print_csc_full_matrix ( sparse_csc_t *A );
int sparse_print_triplet_full_matrix ( sparse_triplet_t *A_triplet );


// ----------------------------------------------------------------------
// Debug with MATLAB
// ----------------------------------------------------------------------
extern int g_matrix_print_format;
#define MATRIX_PRINT_FORMAT_PLAIN 0
#define MATRIX_PRINT_FORMAT_MATLAB 1

// ----------------------------------------------------------------------
// Debug Flags
// ----------------------------------------------------------------------
extern int g_debug_sparse_lu_decomposition;

#endif

