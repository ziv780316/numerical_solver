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

/* DTRSV
	solves one of the systems of equations
	A*x = b, or A**T*x = b,
	where b and x are n element vectors and A is an n by n unit,
	or non-unit, upper or lower triangular matrix.
	No test for singularity or near-singularity is included in this routine. 
	Such tests must be performed before calling this routine.
*/
void dtrsv_( char *triangulr_type, char *trans, char *is_unit_triangular, int *n, double *A, int *lda, double *x, int *incx );

int dense_triangular_solve ( int n, double *A, double *x, bool is_lower_triangular, bool transpose, bool is_unit_triangular )
{
	char type = is_lower_triangular ? 'L' : 'U'; // if is L, will ignore strictly upper triangular part, otherwise
	char tran = transpose ? 'T' : 'N';
	char is_unit = is_unit_triangular ? 'U' : 'N'; // if is U, will regard diagonal element as 1
	int lda = n;
	int incx = 1;

	dtrsv_( &type, &tran, &is_unit, &n, A, &lda, x, &incx );

	return true;
}

/* DSWAP
	interchanges two vectors.
	uses unrolled loops for increments equal to 1.
*/
void dswap_( int *n, double *x, int *incx, double *y, int *incy );

int dense_swap_vector ( int n, double *x, double *y )
{
	int incx = 1;
	int incy = 1;

	dswap_( &n, x, &incx, y, &incy );
	
	return true;
}

/* DGETRS 
	solves a system of linear equations
	A * X = B  or  A**T * X = B
	with a general N-by-N matrix A using the LU factorization computed by DGETRF
*/

int dense_print_vector ( int n, double *x )
{
	for ( int i = 0; i < n; ++i )
	{
		printf( "%.10e\n", x[i] );
	}

	return true;
}


int dense_print_matrix ( int m, int n, double *A )
{
	for ( int i = 0; i < m; ++i )
	{
		for ( int j = 0; j < n; ++j )
		{
			printf( "%.10e ", *(A + j*m + i) );
		}
		printf( "\n" );
	}

	return true;
}

