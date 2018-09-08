#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "matrix_solver.h"

// ============== BLAS and LAPACK Fortran API (column-major array) ============== 
// real matrix
// EX: A[4] = { 1, 2, 3, 4 } 
// A = 1 3
//     2 4
//
// complex matrix
// EX: A[4*2] = { 1, 2, 3, 4, 5, 6, 7, 8 } 
// A = 1+2i 5+6i
//     3+4i 7+8i

/* ZSSCAL
	scales a vector by a complex number
*/
void zscal_( int *n, double *alpha, double *x, int *incx );

// x := alpha*x
int dense_vector_scale ( int n, double *x, double *alpha, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		double r = *alpha;
		for ( int i = 0; i < n; ++i )
		{
			x[i] *= r;
		}
	}
	else
	{
		int incx = 1;
		zscal_( &n, alpha, x, &incx );
	}

	return true;
}

// aij = aii * alpha
int dense_matrix_scale ( int m, int n, double *A, double *alpha, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		double r = *alpha;
		for ( int j = 0; j < n; ++j )
		{
			for ( int i = 0; i < m; ++i )
			{
				*(A + j*m + i) *= r; 
			}
		}
	}
	else
	{
		int incx = 1;
		int col_incr = 2 * m;
		for ( int j = 0; j < n; ++j )
		{
			zscal_( &n, alpha, (A + j*col_incr), &incx );
		}
	}
	return true;
}

// obtain D of A=L+D+U 
int dense_matrix_get_diagonal ( int n, double *A, double *D, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		for ( int i = 0; i < n; ++i )
		{
			D[i] = *(A + i*n + i); 
		}
	}
	else
	{
		int idx;
		for ( int i = 0; i < n; ++i )
		{
			idx = 2 * i;
			D[idx]     = *(A + n * idx + idx); 
			D[idx + 1] = *(A + n * idx + idx + 1); 
		}
	}
	return true;
}

// main diagonal addition, aii = aii + alpha
int dense_diagonal_addition ( int n, double *A, double *alpha, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		double r = *alpha;
		for ( int i = 0; i < n; ++i )
		{
			*(A + i*n + i) += r; 
		}
	}
	else
	{
		int col_incr = 2 * n;
		double real = alpha[0];
		double imag = alpha[1];
		for ( int i = 0; i < n; ++i )
		{
			*(A + i*col_incr + (2 * i)) 	+= real; // real
			*(A + i*col_incr + (2 * i) + 1) += imag; // imag
		}
	}
	return true;
}

/* ZDOTU, ZDOTC
	forms the dot product of two complex vectors (result is complex*16 number, 2 double)
	ZDOTU = X^T * Y
	ZDOTC = X^H * Y
*/
complex_t zdotu_( int *n, double *x, int *incx, double *y, int *incy );
complex_t zdotc_( int *n, double *x, int *incx, double *y, int *incy );

// val = x . y
int dense_vector_inner_product ( int n, double *x, double *y, double *val, bool conjugate, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		double sum = 0.0;
		for ( int i = 0; i < n; ++i )
		{
			sum += x[i] * y[i];
		}
		*val = sum;
	}
	else
	{
		int incx = 1;
		int incy = 1;
		complex_t dot_result;
		if ( conjugate )
		{
			dot_result = zdotc_( &n, x, &incx, y, &incy );
		}
		else
		{
			dot_result = zdotu_( &n, x, &incx, y, &incy );
		}
		val[0] = dot_result.real;
		val[1] = dot_result.imag;
	}
	return true;
}

// A = x . yT
int dense_vector_outer_product ( int m, int n, double *x, double *y, double *A, bool conjugate, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		for ( int j = 0; j < n; ++j )
		{
			for ( int i = 0; i < m; ++i )
			{
				*(A + j*m + i) = x[i] * y[j];
			}
		}
	}
	else
	{
		memset( A, 0, sizeof(double) * 2 * m * n );
		complex_t alpha = { .real = 1.0, .imag = 0.0 };
		dense_maxtrix_rank_1_update ( m, n, A, (double *)&alpha, x, y, conjugate, type );
	}

	return true;
}

// vector norm |x|p = (sum |xi|^p)^(1.0/p)
int dense_vector_norm ( int p_norm, int n, double *x, double *val, number_type type )
{
	switch ( p_norm )
	{
		case -1:
		{
			double max = 0;
			for ( int i = 0; i < n; ++i )
			{
				if ( x[i] > max )
				{
					max = x[i];
				}
			}
			*val = max;
			break;
		}

		case 1:
		{
			double sum = 0.0;
			for ( int i = 0; i < n; ++i )
			{
				sum += x[i];
			}
			*val = sum;
			break;
		}

		case 2:
		{
			double sum = 0.0;
			for ( int i = 0; i < n; ++i )
			{
				sum += x[i] * x[i];
			}
			sum = sqrt( sum );
			break;
		}

		default:
		{
			double sum = 0.0;
			for ( int i = 0; i < n; ++i )
			{
				sum += pow( x[i], p_norm ); 
			}
			if ( 0.0 == sum )
			{
				*val = 0;
			}
			else
			{
				*val = exp( (1.0 / p_norm) * log(sum) );
			}

			break;
		}
	}
	return true;
}

/* ZGERU
	performs the rank 1 operation

	A := alpha*x*y**T + A,

	where alpha is a complex scalar, x is an m element vector, y is an n element
	vector and A is an m by n matrix.
*/
void zgeru_( int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *A, int *lda );

// A = A + alpha * (x . yT)
int dense_maxtrix_rank_1_update ( int m, int n, double *A, double *alpha, double *x, double *y, bool conjugate, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		double r = *alpha;
		for ( int j = 0; j < n; ++j )
		{
			for ( int i = 0; i < m; ++i )
			{
				*(A + j*m + i) += r * (x[i] * y[j]);
			}
		}
	}
	else
	{
		int incx = 1;
		int incy = 1;
		int lda  = m;
		if ( conjugate )
		{
			zgerc_( &m, &n, alpha, x, &incx, y, &incy, A, &lda );
		}
		else
		{
			zgeru_( &m, &n, alpha, x, &incx, y, &incy, A, &lda );
		}
	}
	return true;
}

/* DGEMV
	performs one of the matrix-vector operations
	y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y,
	where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
*/
void dgemv_( char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *x, int *incx, double *beta, double *y, int *incy );

int dense_matrix_vector_multiply ( int m, int n, double *alpha, double *A, double *x, double *beta, double *y, bool transpose, number_type type )
{
	char tran = transpose ? 'T' : 'N';
	int lda = m;
	int incx = 1;
	int incy = 1;

	dgemv_( &tran, &m, &n, alpha, A, &lda, x, &incx, beta, y, &incy );

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

int dense_matrix_matrix_multiply ( int m, int n, int k, double *alpha, double *A, double *B, double *beta, double *C, bool a_transpose, bool b_transpose, number_type type )
{
	char tran_a = a_transpose ? 'T' : 'N';
	char tran_b = b_transpose ? 'T' : 'N';
	int lda = a_transpose ? k : m; 
	int ldb = b_transpose ? n : k; 
	int ldc = m; 

	dgemm_( &tran_a, &tran_b, &m, &n, &k, alpha, A, &lda, B, &ldb, beta, C, &ldc );

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

int dense_triangular_solve ( int n, double *A, double *x, bool is_lower_triangular, bool transpose, bool is_unit_triangular, number_type type )
{
	char triangle_type = is_lower_triangular ? 'L' : 'U'; // if is L, will ignore strictly upper triangular part, otherwise
	char tran = transpose ? 'T' : 'N';
	char is_unit = is_unit_triangular ? 'U' : 'N'; // if is U, will regard diagonal element as 1
	int lda = n;
	int incx = 1;

	dtrsv_( &triangle_type, &tran, &is_unit, &n, A, &lda, x, &incx );

	return true;
}

/* DSWAP
	interchanges two vectors.
	uses unrolled loops for increments equal to 1.
*/
void dswap_( int *n, double *x, int *incx, double *y, int *incy );

int dense_swap_vector ( int n, double *x, double *y, number_type type )
{
	int incx = 1;
	int incy = 1;

	dswap_( &n, x, &incx, y, &incy );
	
	return true;
}

/* DGETRF DGETRF computes an LU factorization of a general M-by-N matrix A 
	using partial pivoting with row interchanges.
	The factorization has the form
	A = P * L * U
	where P is a permutation matrix, 
	L is lower triangular with unit diagonal elements (lower trapezoidal if m > n),
	and U is upper triangular (upper trapezoidal if m < n).
*/
void dgetrf_ ( int *m, int *n, double *A, int *lda, int*ipiv, int *info );

int dense_lu_factor ( int n, double *A, int *p, number_type type )
{
	int lda = n;
	int info;
	dgetrf_( &n, &n, A, &lda, p, &info );

	return (0 == info); // info = 0 means success
}

/* DGETRS 
	solves a system of linear equations
	A * X = B  or  A**T * X = B
	with a general N-by-N matrix A using the LU factorization computed by DGETRF
*/
void dgetrs_( char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *x, int *ldb, int *info );

int dense_solve ( int n, double *A, double *x, int *ipiv, bool transpose, number_type type )
{
	char tran = transpose ? 'T' : 'N';
	int nrhs = 1;
	int lda = n;
	int ldb = n;
	int info;

	dgetrs_( &tran, &n, &nrhs, A, &lda, ipiv, x, &ldb, &info );
	
	return (0 == info); // info = 0 means success
}

int dense_factor_and_solve ( int n, double *A, double *x, bool transpose, number_type type )
{
	int *p = (int *) malloc( sizeof(int) * n );
	int success;

	success = dense_lu_factor( n, A, p, type );
	if ( !success )
	{
		return 0;
	}
	success = dense_solve( n, A, x, p, transpose, type );
	if ( !success )
	{
		return 0;
	}

	free(p);

	return true;
}

/* DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
   This method inverts U and then computes inv(A) by solving the system inv(A)*L = inv(U) for inv(A).
*/
void dgetri_ ( int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info );

int dense_matrix_inverse ( int n, double *A, int *p, number_type type )
{
	int lda = n;
	int lwork = -1;
	int info;
	double optima_lwork;
	double *work;
	bool success;

	// need LU factor before invese
	success = dense_lu_factor( n, A, p, type );
	if ( !success )
	{
		return 0;
	}

	// query suitable worksapce size first
	dgetri_( &n, A, &lda, p, &optima_lwork, &lwork, &info );

	// allocate temperal memory for optimizing performance
	lwork = (int) optima_lwork;
	work = (double *) malloc (sizeof(double) * lwork);
	//fprintf( stderr, "[matrix info] %s: n=%d optimized_lwork=%d\n", __func__, n, lwork );

	// inverse A
	dgetri_( &n, A, &lda, p, work, &lwork, &info );

	free( work );
	
	return (0 == info); // info = 0 means success, info=i > 0 means singular and A(i,i) is 0 (singular node)
}

int dense_print_vector ( int n, double *x, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		for ( int i = 0; i < n; ++i )
		{
			printf( "%.10e\n", x[i] );
		}
	}
	else
	{
		for ( int i = 0; i < 2 * n; i += 2 )
		{
			printf( "%.10e + i*%.10e\n", x[i], x[i + 1] );
		}
	}

	return true;
}

int dense_print_vector_i ( int n, int *x, number_type type )
{
	for ( int i = 0; i < n; ++i )
	{
		printf( "%d\n", x[i] );
	}

	return true;
}


int dense_print_matrix ( int m, int n, double *A, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		for ( int i = 0; i < m; ++i )
		{
			for ( int j = 0; j < n; ++j )
			{
				printf( "%.10e ", *(A + j*m + i) );
			}
			printf( "\n" );
		}
	}
	else
	{
		int col_incr = 2 * m;
		for ( int i = 0; i < col_incr; i += 2 )
		{
			for ( int j = 0; j < n; ++j )
			{
				printf( "%.10e+i*%.10e ", *(A + j*col_incr + i), *(A + j*col_incr + i + 1) );
			}
			printf( "\n" );
		}
	}


	return true;
}

int dense_print_matrix_LU ( int n, double *A, number_type type )
{
	printf( "L=\n" ); // assume L is unit lower triangular
	for ( int i = 0; i < n; ++i )
	{
		for ( int j = 0; j <= i; ++j )
		{
			if ( i == j )
			{
				printf( "%.10e", 1.0 );
			}
			else
			{
				printf( "%.10e ", *(A + j*n + i) );
			}
		}
		printf( "\n" );
	}

	printf( "U=\n" );
	for ( int i = 0; i < n; ++i )
	{
		for ( int j = i; j < n; ++j )
		{
			printf( "%.10e ", *(A + j*n + i) );
		}
		printf( "\n" );
	}

	return true;
}
