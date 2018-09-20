#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "matrix_solver.h"

static void iswap ( int *x, int *y )
{
	// XOR swap does not need buffer
	*x = *x ^ *y;
	*y = *x ^ *y;
	*x = *x ^ *y;
}

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
void zscal ( int *n, double *alpha, double *x, int *incx );

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
		zscal( &n, alpha, x, &incx );
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
			zscal( &n, alpha, (A + j*col_incr), &incx );
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
complex_t zdotu ( int *n, double *x, int *incx, double *y, int *incy );
complex_t zdotc ( int *n, double *x, int *incx, double *y, int *incy );

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
			dot_result = zdotc( &n, x, &incx, y, &incy );
		}
		else
		{
			dot_result = zdotu( &n, x, &incx, y, &incy );
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
			if ( REAL_NUMBER == type )
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
			}
			else
			{
				double max = 0.0;
				double z_val;
				for ( int i = 0; i < 2 * n; i += 2 )
				{
					z_val = x[i]*x[i] + x[i+1]*x[i+1];
					if ( z_val > max )
					{
						max = z_val;
					}
				}
				max = sqrt( max );
				*val = max;
			}
			break;
		}

		case 1:
		{
			if ( REAL_NUMBER == type )
			{
				double sum = 0.0;
				for ( int i = 0; i < n; ++i )
				{
					sum += x[i];
				}
				*val = sum;
			}
			else
			{
				double sum = 0.0;
				for ( int i = 0; i < 2 * n; i += 2 )
				{
					sum += sqrt( x[i]*x[i] + x[i+1]*x[i+1] );
				}
				*val = sum;
			}
			break;
		}

		case 2:
		{
			if ( REAL_NUMBER == type )
			{
				double sum = 0.0;
				for ( int i = 0; i < n; ++i )
				{
					sum += x[i] * x[i];
				}
				sum = sqrt( sum );
				*val = sum;
			}
			else
			{
				double sum = 0.0;
				for ( int i = 0; i < 2 * n; i += 2 )
				{
					sum += x[i]*x[i] + x[i+1]*x[i+1];
				}
				sum = sqrt( sum );
				*val = sum;
			}
			break;
		}

		default:
		{
			if ( REAL_NUMBER == type )
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
			}
			else
			{
				double sum = 0.0;
				double z_val;
				for ( int i = 0; i < 2 * n; i += 2 )
				{
					z_val = sqrt(x[i]*x[i] + x[i+1]*x[i+1]);
					sum += pow( z_val, p_norm ); 
				}
				if ( 0.0 == sum )
				{
					*val = 0;
				}
				else
				{
					*val = exp( (1.0 / p_norm) * log(sum) );
				}
			}

			break;
		}
	}
	return true;
}

/* ZGERU, ZGERC
	performs the rank 1 operation

	A := alpha*x*(y**T) + A or A := alpha*x*(y**H) + A

	where alpha is a complex scalar, x is an m element vector, y is an n element
	vector and A is an m by n matrix.
*/
void zgeru ( int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *A, int *lda );
void zgerc ( int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *A, int *lda );

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
			// y conjugate
			zgerc( &m, &n, alpha, x, &incx, y, &incy, A, &lda );
		}
		else
		{
			zgeru( &m, &n, alpha, x, &incx, y, &incy, A, &lda );
		}
	}
	return true;
}

/* DGEMV, ZGEMV
	performs one of the matrix-vector operations
	y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y or y := alpha*A**H*x + beta*y
	where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.

	TRANS = 'N' or 'n' y := alpha*A*x + beta*y.
	TRANS = 'T' or 't' y := alpha*A**T*x + beta*y.
	TRANS = 'C' or 'c' y := alpha*A**H*x + beta*y.
*/
void dgemv ( char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *x, int *incx, double *beta, double *y, int *incy );
void zgemv ( char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *x, int *incx, double *beta, double *y, int *incy );

int dense_matrix_vector_multiply ( int m, int n, double *alpha, double *A, double *x, double *beta, double *y, transpose_type trans_type, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		char tran = 'N';
		int lda = m;
		int incx = 1;
		int incy = 1;

		if ( TRANS_NORMAL == trans_type )
		{
			tran = 'T';
		}
		
		dgemv( &tran, &m, &n, alpha, A, &lda, x, &incx, beta, y, &incy );
	}
	else
	{
		char tran = 'N';
		int lda = m;
		int incx = 1;
		int incy = 1;

		if ( TRANS_NORMAL == trans_type )
		{
			tran = 'T';
		}
		else if ( TRANS_CONJUGATE == trans_type )
		{
			tran = 'C';
		}

		zgemv( &tran, &m, &n, alpha, A, &lda, x, &incx, beta, y, &incy );
	}
	return true;
}

/* DGEMM, ZGEMM
	performs one of the matrix-matrix operations
	C := alpha*op(A)*op(B) + beta*C,
	where  op(X) is one of op(X) = X or op(X) = X**T or op(X) = X**H
	alpha and beta are scalars, and A, B and C are matrices
	op(A) an m by k matrix, op(B) a k by n matrix and C an m by n matrix.

	TRANS = 'N' or 'n' op(A)=A
	TRANS = 'T' or 't' op(A)=A**T
	TRANS = 'C' or 'c' op(A)=A**H
*/
void dgemm ( char *a_transpose, char *b_transpose, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc );
void zgemm ( char *a_transpose, char *b_transpose, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc );

int dense_matrix_matrix_multiply ( int ma, int na, int mb, int nb, double *alpha, double *A, double *B, double *beta, double *C, transpose_type a_transpose, transpose_type b_transpose, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		char tran_a = 'N';
		char tran_b = 'N';
		int lda = ma; 
		int ldb = mb; 
		int ldc;

		if ( TRANS_NORMAL == a_transpose )
		{
			tran_a = 'T';
			iswap( &ma, &na );
		}
		if ( TRANS_NORMAL == b_transpose )
		{
			tran_b = 'T';
			iswap( &mb, &nb );
		}

		if ( na != mb )
		{
			fprintf( stderr, "[Error] matrix dimension mismatch (A=%dx%d B=%dx%d)\n", ma, na, mb, nb );
			abort();
		}

		ldc = ma;

		dgemm( &tran_a, &tran_b, &ma, &nb, &na, alpha, A, &lda, B, &ldb, beta, C, &ldc );
	}
	else
	{
		char tran_a = 'N';
		char tran_b = 'N';
		int lda = ma; 
		int ldb = mb; 
		int ldc;

		if ( TRANS_NORMAL == a_transpose )
		{
			tran_a = 'T';
			iswap( &ma, &na );
		}
		if ( TRANS_NORMAL == b_transpose )
		{
			tran_b = 'T';
			iswap( &mb, &nb );
		}
		if ( TRANS_CONJUGATE == a_transpose )
		{
			tran_a = 'C';
			iswap( &ma, &na );
		}
		if ( TRANS_CONJUGATE == b_transpose )
		{
			tran_b = 'C';
			iswap( &mb, &nb );
		}

		ldc = ma;

		zgemm( &tran_a, &tran_b, &ma, &nb, &na, alpha, A, &lda, B, &ldb, beta, C, &ldc );
	}
	return true;
}

/* DTRSV, ZTRSV
	solves one of the systems of equations
	A*x = b, or A**T*x = b,
	where b and x are n element vectors and A is an n by n unit,
	or non-unit, upper or lower triangular matrix.
	No test for singularity or near-singularity is included in this routine. 
	Such tests must be performed before calling this routine.
*/
void dtrsv ( char *trig, char *tran, char *is_unit_triangular, int *n, double *A, int *lda, double *x, int *incx );
void ztrsv ( char *trig, char *tran, char *is_unit_triangular, int *n, double *A, int *lda, double *x, int *incx );

int dense_triangular_solve ( int n, double *A, double *x, triangular_type triangle_type, transpose_type transpose, number_type type )
{
	char trig = ((TRIG_LOWER == triangle_type) ||
		     (TRIG_LOWER_UNIT == triangle_type)) ? 'L' : 'U'; // if is L, will ignore strictly upper triangular part, otherwise
	char tran = 'N';
	char is_unit = ((TRIG_LOWER_UNIT == triangle_type) || 
			(TRIG_UPPER_UNIT == triangle_type)) ? 'U' : 'N'; // if is U, will regard diagonal element as 1
	int lda = n;
	int incx = 1;

	if ( TRANS_NORMAL == transpose )
	{
		tran = 'T';
	}
	else if ( TRANS_CONJUGATE == transpose )
	{
		tran = 'C';
	}

	if ( type == REAL_NUMBER )
	{
		dtrsv( &trig, &tran, &is_unit, &n, A, &lda, x, &incx );
	}
	else
	{
		ztrsv( &trig, &tran, &is_unit, &n, A, &lda, x, &incx );
	}
	return true;
}

/* DSWAP, ZSWAP
	interchanges two vectors.
	uses unrolled loops for increments equal to 1.
*/
void dswap ( int *n, double *x, int *incx, double *y, int *incy );
void zswap ( int *n, double *x, int *incx, double *y, int *incy );

int dense_swap_vector ( int n, double *x, double *y, number_type type )
{
	int incx = 1;
	int incy = 1;

	if ( REAL_NUMBER == type )
	{
		dswap( &n, x, &incx, y, &incy );
	}
	else
	{
		zswap( &n, x, &incx, y, &incy );
	}
	
	return true;
}

/* DGETRF, ZGETRF 
	computes an LU factorization of a general M-by-N matrix A 
	using partial pivoting with row interchanges.
	The factorization has the form
	A = P * L * U
	where P is a permutation matrix, 
	L is lower triangular with unit diagonal elements (lower trapezoidal if m > n),
	and U is upper triangular (upper trapezoidal if m < n).
*/
void dgetrf ( int *m, int *n, double *A, int *lda, int *ipiv, int *info );
void zgetrf ( int *m, int *n, double *A, int *lda, int *ipiv, int *info );

int dense_lu_factor ( int n, double *A, int *p, number_type type )
{
	int lda = n;
	int info;

	if ( REAL_NUMBER == type )
	{
		dgetrf( &n, &n, A, &lda, p, &info );
	}
	else
	{
		zgetrf( &n, &n, A, &lda, p, &info );
	}

	return (0 == info); // info = 0 means success
}

/* DGETRS, ZGETRS
	solves a system of linear equations
	A * X = B  or  A**T * X = B or A**H X = B
	with a general N-by-N matrix A using the LU factorization computed by DGETRF
*/
void dgetrs ( char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *x, int *ldb, int *info );
void zgetrs ( char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *x, int *ldb, int *info );

int dense_solve ( int n, int nrhs, double *A, double *x, int *ipiv, transpose_type transpose, number_type type )
{
	char tran = 'N';
	int lda = n;
	int ldb = n;
	int info;

	if ( TRANS_NORMAL == transpose )
	{
		tran = 'T';
	}
	else if ( TRANS_CONJUGATE == transpose )
	{
		tran = 'C';
	}

	if ( REAL_NUMBER == type )
	{
		dgetrs( &tran, &n, &nrhs, A, &lda, ipiv, x, &ldb, &info );
	}
	else
	{
		zgetrs( &tran, &n, &nrhs, A, &lda, ipiv, x, &ldb, &info );
	}
	
	return (0 == info); // info = 0 means success
}

int dense_factor_and_solve ( int n, int nrhs, double *A, double *x, transpose_type tran, number_type type )
{
	int *p = (int *) malloc( sizeof(int) * n );
	int success;

	success = dense_lu_factor( n, A, p, type );
	if ( !success )
	{
		return 0;
	}
	success = dense_solve( n, nrhs, A, x, p, tran, type );
	if ( !success )
	{
		return 0;
	}

	free(p);

	return true;
}

/* DGETRI ZGETRI
	computes the inverse of a matrix using the LU factorization computed by DGETRF.
   	This method inverts U and then computes inv(A) by solving the system inv(A)*L = inv(U) for inv(A).
*/
void dgetri ( int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info );
void zgetri ( int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info );

int dense_matrix_inverse ( int n, double *A, int *p, number_type type )
{
	int lda = n;
	int lwork = -1;
	int info;
	double optima_lwork[2]; // need at least 2 double for zgetri to prevent stack smash (there is complex16 work(1) access in zgetri)
	double *work;
	bool success;

	// need LU factor before invese
	success = dense_lu_factor( n, A, p, type );
	if ( !success )
	{
		return 0;
	}

	// query suitable worksapce size first
	if ( REAL_NUMBER == type )
	{
		dgetri( &n, A, &lda, p, optima_lwork, &lwork, &info );
	}
	else
	{
		zgetri( &n, A, &lda, p, optima_lwork, &lwork, &info );
	}
	//fprintf( stderr, "[matrix info] %s: n=%d optimized_lwork=%d\n", __func__, n, (int)optima_lwork[0] );

	// inverse A
	lwork = (int) optima_lwork[0];
	if ( REAL_NUMBER == type )
	{
		// allocate temperal memory for optimizing performance
		work = (double *) malloc (sizeof(double) * lwork);
		dgetri( &n, A, &lda, p, work, &lwork, &info );
	}
	else
	{
		work = (double *) malloc (sizeof(double) * 2 * lwork);
		zgetri( &n, A, &lda, p, work, &lwork, &info );
	}

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
	if ( REAL_NUMBER == type )
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
	}
	else
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
					printf( "%.10e+i%.10e ", *(A + j*2*n + 2*i), *(A + j*2*n + 2*i + 1) );
				}
			}
			printf( "\n" );
		}

		printf( "U=\n" );
		for ( int i = 0; i < n; ++i )
		{
			for ( int j = i; j < n; ++j )
			{
				printf( "%.10e+i%.10e ", *(A + j*2*n + 2*i), *(A + j*2*n + 2*i + 1) );
			}
			printf( "\n" );
		}
	}

	return true;
}
