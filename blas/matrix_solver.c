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
	char tran = transpose ? 'Y' : 'N';
	int lda = m; 
	int incx = 1;
	int incy = 1;

	dgemv_( &tran, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );

	return true;
}
