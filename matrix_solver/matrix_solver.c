#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "matrix_solver.h"

int g_matrix_print_format = MATRIX_PRINT_FORMAT_PLAIN;

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
					if ( fabs(x[i]) > max )
					{
						max = fabs(x[i]);
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
					sum += fabs(x[i]);
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
					sum += pow( fabs(x[i]), p_norm ); 
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
			exit(1);
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

/* 1. DGETRF, ZGETRF 
	computes an LU factorization of a general M-by-N matrix A 
	using partial pivoting with row interchanges.
	The factorization has the form
	A = P * L * U
	where P is a permutation matrix, 
	L is lower triangular with unit diagonal elements (lower trapezoidal if m > n),
	and U is upper triangular (upper trapezoidal if m < n).
	ipiv is P**T

   2. DOPTRF, ZPOTRF
	computes the Cholesky factorization of a real symmetric positive definite matrix A.
	The factorization has the form
	A = U**T * U,  if UPLO = 'U', or
	A = L  * L**T,  if UPLO = 'L',
	where U is an upper triangular matrix and L is lower triangular.

   3. DSYTRF, ZSYTRF
	computes the factorization of a real symmetric matrix A using
	the Bunch-Kaufman diagonal pivoting method.  The form of the factorization is
	A = U*D*U**T  or  A = L*D*L**T
	where U (or L) is a product of permutation and unit upper (lower)
	triangular matrices, and D is symmetric and block diagonal with
	1-by-1 and 2-by-2 diagonal blocks.
*/

void dgetrf ( int *m, int *n, double *A, int *lda, int *ipiv, int *info );
void zgetrf ( int *m, int *n, double *A, int *lda, int *ipiv, int *info );
void dpotrf ( char *uplo, int *n, double *A, int *lda, int *info );
void zpotrf ( char *uplo, int *n, double *A, int *lda, int *info );
void dsytrf ( char *uplo, int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info );
void zsytrf ( char *uplo, int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info );

int dense_lu_factor ( int n, double *A, int *pinv, factorization_type factor_method, number_type type )
{
	int lda = n;
	int info;

	if ( FACTOR_LU_RIGHT_LOOKING == factor_method )
	{
		if ( REAL_NUMBER == type )
		{
			dgetrf( &n, &n, A, &lda, pinv, &info );
		}
		else
		{
			zgetrf( &n, &n, A, &lda, pinv, &info );
		}
	}
	else if ( FACTOR_LU_CHOLESKY == factor_method )
	{
		char uplo = 'L';
		if ( REAL_NUMBER == type )
		{
			dpotrf( &uplo, &n, A, &lda, &info );
		}
		else
		{
			zpotrf( &uplo, &n, A, &lda, &info );
		}
	}
	else if ( FACTOR_LU_BUNCH_KAUFMAN == factor_method )
	{
		char uplo = 'L';
		int lwork = -1;
		double optima_lwork[2]; // need at least 2 double for zgetri to prevent stack smash (there is complex16 work(1) access in zsytrf)
		double *work;

		// query lwork
		if ( REAL_NUMBER == type )
		{
			dsytrf ( &uplo, &n, A, &lda, pinv, optima_lwork, &lwork, &info );
		}
		else
		{
			zsytrf ( &uplo, &n, A, &lda, pinv, optima_lwork, &lwork, &info );
		}

		//fprintf( stderr, "[matrix info] %s: n=%d optima_lwork=%d\n", __func__, n, (int)optima_lwork[0] );

		// factorization
		lwork = (int) optima_lwork[0];
		if ( REAL_NUMBER == type )
		{
			// allocate temperal memory for optimizing performance
			work = (double *) malloc (sizeof(double) * lwork);
			dsytrf ( &uplo, &n, A, &lda, pinv, work, &lwork, &info );
		}
		else
		{
			// allocate temperal memory for optimizing performance
			work = (double *) malloc (sizeof(double) * 2 * lwork);
			zsytrf ( &uplo, &n, A, &lda, pinv, work, &lwork, &info );
		}

		free( work );
	}

	return (0 == info); // info = 0 means success
}

/* 1. DGETRS, ZGETRS
	solves a system of linear equations
	A * X = B  or  A**T * X = B or A**H X = B
	with a general N-by-N matrix A using the LU factorization computed by DGETRF
	
   2. DPOTRS, ZPOTRS
	solves a system of linear equations A*X = B with a symmetric
	positive definite matrix A using the Cholesky factorization
	A = U**T*U or A = L*L**T computed by DPOTRF.

   3. DSYTRS, ZSYTRS
	solves a system of linear equations A*X = B with a real
	symmetric matrix A using the factorization A = U*D*U**T or
	A = L*D*L**T computed by DSYTRF.

*/
void dgetrs ( char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *x, int *ldb, int *info );
void zgetrs ( char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *x, int *ldb, int *info );
void dpotrs ( char *uplo, int *n, int *nrhs, double *A, int *lda, double *x, int *ldb, int *info );
void zpotrs ( char *uplo, int *n, int *nrhs, double *A, int *lda, double *x, int *ldb, int *info );
void dsytrs ( char *uplo, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *x, int *ldb, int *info );
void zsytrs ( char *uplo, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *x, int *ldb, int *info );

int dense_solve ( int n, int nrhs, double *A, double *x, int *ipiv, factorization_type factor_method, transpose_type transpose, number_type type )
{
	int lda = n;
	int ldb = n;
	int info;

	if ( FACTOR_LU_RIGHT_LOOKING == factor_method )
	{
		char tran = 'N';

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
	}
	else if ( FACTOR_LU_CHOLESKY == factor_method )
	{
		char uplo = 'L';

		if ( REAL_NUMBER == type )
		{
			dpotrs ( &uplo, &n, &nrhs, A, &lda, x, &ldb, &info );
		}
		else
		{
			zpotrs ( &uplo, &n, &nrhs, A, &lda, x, &ldb, &info );
		}
	}
	else if ( FACTOR_LU_BUNCH_KAUFMAN == factor_method )
	{
		char uplo = 'L';

		if ( REAL_NUMBER == type )
		{
			dsytrs( &uplo, &n, &nrhs, A, &lda, ipiv, x, &ldb, &info );
		}
		else
		{
			zsytrs( &uplo, &n, &nrhs, A, &lda, ipiv, x, &ldb, &info );
		}
	}

	return (0 == info); // info = 0 means success
}

// eval |A|
double dense_eval_det ( int n, double *A, int *pinv, factorization_type factor_method, number_type type )
{
	double det_a = NAN;
	if ( REAL_NUMBER == type )
	{
		if ( FACTOR_LU_BUNCH_KAUFMAN == factor_method ) // LDL' for all symmetry matrix
		{
		}
		else if ( FACTOR_LU_CHOLESKY == factor_method ) // LL' for postive definite symmetric matrix
		{
		}
		else
		{
			// A = P * L * U, |A| = |P| * |L| * |U| = |P| * |U| = -|P'| * |U|

			// |P'|
			int swap_cnt = 0;
			double det_p;
			double det_pinv;
			for ( int i = 1; i <= n; ++i )
			{
				if ( i != pinv[i - 1] ) 
				{
					++swap_cnt;
				}
			}
			if ( 1 == (swap_cnt % 2) )
			{
				det_pinv = -1;
				det_p = -1;
			}
			else
			{
				det_pinv = 1;
				det_p = 1;
			}

			// |U|
			double det_u = 1.0;
			for ( int i = 0; i < n; ++i )
			{
				det_u *= *(A + i*n + i);
			}

			// |A| = |P| * |U|
			det_a = det_p * det_u;
		}
	}
	else
	{
		fprintf( stderr, "[Error] cannot support evaluate complex matrix |A|\n" );
		exit(1);
	}

	return det_a;
}

int dense_factor_and_solve ( int n, int nrhs, double *A, double *x, int *pinv, factorization_type factor_method, transpose_type tran, number_type type )
{
	int success;

	success = dense_lu_factor( n, A, pinv, factor_method, type );
	if ( !success )
	{
		return 0;
	}
	success = dense_solve( n, nrhs, A, x, pinv, factor_method, tran, type );
	if ( !success )
	{
		return 0;
	}

	return true;
}

/* DGETRI ZGETRI
	computes the inverse of a matrix using the LU factorization computed by DGETRF.
   	This method inverts U and then computes inv(A) by solving the system inv(A)*L = inv(U) for inv(A).
*/
void dgetri ( int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info );
void zgetri ( int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info );
void dpotri ( char *uplo, int *n, double *A, int *lda, int *info );
void zpotri ( char *uplo, int *n, double *A, int *lda, int *info );
void dsytri ( char *uplo, int *n, double *A, int *lda, int *ipiv, double *work, int *info );
void zsytri ( char *uplo, int *n, double *A, int *lda, int *ipiv, double *work, int *info );

int dense_matrix_inverse ( int n, double *A, int *pinv, factorization_type factor_method, number_type type )
{
	int lda = n;
	int lwork = -1;
	int info;
	double optima_lwork[2]; // need at least 2 double for zgetri to prevent stack smash (there is complex16 work(1) access in zgetri)
	double *work;
	bool success;

	// need factorization before invese
	success = dense_lu_factor( n, A, pinv, factor_method, type );
	if ( !success )
	{
		return 0;
	}

	if ( FACTOR_LU_RIGHT_LOOKING == factor_method )
	{
		// query suitable worksapce size first, lwork = -1 means query optima workspace size (store in optima_lwork[0])
		if ( REAL_NUMBER == type )
		{
			dgetri( &n, A, &lda, pinv, optima_lwork, &lwork, &info );
		}
		else
		{
			zgetri( &n, A, &lda, pinv, optima_lwork, &lwork, &info );
		}
		//fprintf( stderr, "[matrix info] %s: n=%d optima_lwork=%d\n", __func__, n, (int)optima_lwork[0] );

		// inverse A
		lwork = (int) optima_lwork[0];
		if ( REAL_NUMBER == type )
		{
			// allocate temperal memory for optimizing performance
			work = (double *) malloc (sizeof(double) * lwork);
			dgetri( &n, A, &lda, pinv, work, &lwork, &info );
		}
		else
		{
			work = (double *) malloc (sizeof(double) * 2 * lwork);
			zgetri( &n, A, &lda, pinv, work, &lwork, &info );
		}

		free( work );
	}
	else if ( FACTOR_LU_CHOLESKY == factor_method )
	{
		// query suitable worksapce size first, lwork = -1 means query optima workspace size (store in optima_lwork[0])
		char uplo = 'L';
		if ( REAL_NUMBER == type )
		{
			dpotri( &uplo, &n, A, &lda, &info );
		}
		else
		{
			zpotri( &uplo, &n, A, &lda, &info );
		}

		// only inverse lower triangle, so add upper manual
		for ( int j = 0; j < n; ++j )
		{
			for ( int i = j + 1; i < n; ++i )
			{
				*(A + i*n + j) = *(A + j*n + i);
			}
		}
	}
	else if ( FACTOR_LU_BUNCH_KAUFMAN == factor_method )
	{
		// query suitable worksapce size first, lwork = -1 means query optima workspace size (store in optima_lwork[0])
		char uplo = 'L';
		if ( REAL_NUMBER == type )
		{
			work = (double *) malloc (sizeof(double) * n );
			dsytri( &uplo, &n, A, &lda, pinv, work, &info );
		}
		else
		{
			work = (double *) malloc (sizeof(double) * 2 * n );
			zsytri( &uplo, &n, A, &lda, pinv, work, &info );
		}

		free( work );

		// only inverse lower triangle, so add upper manual
		for ( int j = 0; j < n; ++j )
		{
			for ( int i = j + 1; i < n; ++i )
			{
				*(A + i*n + j) = *(A + j*n + i);
			}
		}
	}
	
	return (0 == info); // info = 0 means success, info=i > 0 means singular and A(i,i) is 0 (singular node)
}

int dense_matrix_norm ( int p_norm, int n, int m, double *A, double *val, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		switch ( p_norm )
		{
			// max
			case -1:
				{
					double norm = -DBL_MAX;
					double sum;
					for ( int i = 0; i < m; ++i )
					{
						sum = 0.0;
						for ( int j = 0; j < n; ++j )
						{
							sum += fabs(*(A + j*m + i));
						}
						norm = fmax( norm, fabs(sum) );
					}
					*val = norm;
					break;
				}

			case 1:
			case 2:
			default:
				fprintf( stderr, "[Error] cannot support matrix %d norm\n", p_norm );
				exit(1);
				break;
		}
	}
	else
	{
		fprintf( stderr, "[Error] cannot support complex matrix norm\n" );
		exit(1);
	}


	return true;
}

int dense_print_vector ( int n, double *x, number_type type )
{
	if ( MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format )
	{
		printf( "[...\n" );
	}
	if ( REAL_NUMBER == type )
	{
		for ( int i = 0; i < n; ++i )
		{
			printf( "%.10e%s\n", x[i], ((MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format) ? ";..." : "") );
		}
	}
	else
	{
		for ( int i = 0; i < 2 * n; i += 2 )
		{
			printf( "%.10e + i*%.10e%s\n", x[i], x[i + 1], ((MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format) ? ";..." : "") );
		}
	}
	if ( MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format )
	{
		printf( "];\n" );
	}

	return true;
}

int dense_print_vector_i ( int n, int *x )
{
	if ( MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format )
	{
		printf( "[..." );
	}
	for ( int i = 0; i < n; ++i )
	{
		printf( "%d%s\n", x[i], ((MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format) ? ";..." : "") );
	}
	if ( MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format )
	{
		printf( "];\n" );
	}

	return true;
}


int dense_print_matrix ( int m, int n, double *A, number_type type )
{
	if ( MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format )
	{
		printf( "[...\n" );
	}
	if ( REAL_NUMBER == type )
	{
		for ( int i = 0; i < m; ++i )
		{
			for ( int j = 0; j < n; ++j )
			{
				printf( "%.10e ", *(A + j*m + i) );
			}
			printf( "%s\n", ((MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format) ? ";..." : "") );
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
			printf( "%s\n", ((MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format) ? ";..." : "") );
		}
	}
	if ( MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format )
	{
		printf( "];\n" );
	}

	return true;
}

int dense_print_matrix_perm ( int n, int *p )
{
	if ( MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format )
	{
		printf( "[...\n" );
	}

	int *idx = (int *) malloc ( sizeof(int) * (n + 1) );
	for ( int i = 1; i <= n; ++i )
	{
		idx[i] = i;
	}

	for ( int i = 1; i <= n; ++i )
	{
		if ( i != p[i - 1] ) 
		{
			iswap( &(idx[i]), &(idx[ p[i - 1] ]) );
		}
	}
	
	//printf( "permutation=\n");
	//for ( int i = 1; i <= n; ++i )
	//{
	//	printf( "%d ", idx[i] );
	//}
	//printf("\n");

	for ( int i = 1; i <= n; ++i )
	{
		for ( int j = 1; j <= n; ++j )
		{
			if ( j == idx[i] )
			{
				printf( "1 ");
			}
			else
			{
				printf( "0 ");
			}
		}
		printf( "%s\n", ((MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format) ? ";..." : "") );
	}

	if ( MATRIX_PRINT_FORMAT_MATLAB == g_matrix_print_format )
	{
		printf( "];\n" );
	}

	free( idx );
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

int dense_print_matrix_trig ( int n, double *A, triangular_type trig, number_type type )
{
	if ( REAL_NUMBER == type )
	{
		if ( TRIG_LOWER == trig )
		{
			printf( "L=\n" ); 
			for ( int i = 0; i < n; ++i )
			{
				for ( int j = 0; j <= i; ++j )
				{
					printf( "%.10e ", *(A + j*n + i) );
				}
				printf( "\n" );
			}
		}
		else if ( TRIG_LOWER_UNIT == trig )
		{
			printf( "L=\n" ); 
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
		}
		else if ( TRIG_UPPER == trig )
		{
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
		else if ( TRIG_UPPER_UNIT == trig )
		{
			printf( "U=\n" ); 
			for ( int i = 0; i < n; ++i )
			{
				for ( int j = i; j <= n; ++j )
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
		}
	}
	else
	{
		if ( TRIG_LOWER == trig )
		{
			printf( "L=\n" ); 
			for ( int i = 0; i < n; ++i )
			{
				for ( int j = 0; j <= i; ++j )
				{
					printf( "%.10e+i%.10e ", *(A + j*2*n + 2*i), *(A + j*2*n + 2*i + 1) );
				}
				printf( "\n" );
			}
		}
		else if ( TRIG_LOWER_UNIT == trig )
		{
			printf( "L=\n" ); 
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
		}
		else if ( TRIG_UPPER == trig )
		{
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
		else if ( TRIG_UPPER_UNIT == trig )
		{
			printf( "U=\n" );
			for ( int i = 0; i < n; ++i )
			{
				for ( int j = i; j < n; ++j )
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
		}
	}
	return true;
}
