#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

// ============== BLAS Fortran API (column-major array) ============== 
/* DGEMV
	performs one of the matrix-vector operations
	y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y,
	where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
*/
void dgemv_( char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *x, int *incx, double *beta, double *y, int *incy );

int dense_matrix_vector_multiply ( int m, int n, double alpha, double *A, double *x, double beta, double *y, bool transpose )
{
	char tran = transpose ? 'T' : 'N';
	int lda = m;
	int incx = 1;
	int incy = 1;

	dgemv_( &tran, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );

	return true;
}

/* DGEMM  
	performs one of the matrix-matrix operations
	C := alpha*op(A)*op(B) + beta*C,
	where  op(X  is one of op(X) = X or op(X) = X**T,
	alpha and beta are scalars, and A, B and C are matrices
	op(A) an m by k matrix, op(B) a k by n matrix and C an m by n matrix.
*/
void dgemm_( char *a_transpose, char *b_transpose, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc );

int dense_matrix_matrix_multiply ( int m, int n, int k, double alpha, double *A, double *B, double beta, double *C, bool a_transpose, bool b_transpose )
{
	char tran_a = a_transpose ? 'T' : 'N';
	char tran_b = b_transpose ? 'T' : 'N';
	int lda = a_transpose ? k : m; 
	int ldb = b_transpose ? n : k; 
	int ldc = m; 

	dgemm_( &tran_a, &tran_b, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );

	return true;
}
