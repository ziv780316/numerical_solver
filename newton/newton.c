#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <errno.h>

#include "newton.h"
#include "matrix_solver.h"

static void check_user_define_args ( double *x0,
				     void (load_f) (double *x, double*f),
				     void (load_jacobian) (double *x, double*J) );
static void newton_initialize ( int n, double *x, double *x0, bool random_initial );
static void broyden_update ( int n, double *J, double *df, double *dx, bool debug );
static void broyden_update_sherman_morrison ( int n, double *J, double *df, double *dx, bool debug );
static bool __bypass_check ( int n, double *x, double *f, double *dx, double bypass_rtol, double bypass_atol );

bool newton_solve ( newton_iterative_type iterative_type, 
		    newton_damped_type damped_type,
		    newton_rescue_type rescue_type,
		    newton_derivative_type diff_type,
		    int n,
		    double *x0,
		    double *x_result,
		    double *f_result,
		    void (load_f) (double *x, double*f),
		    void (load_jacobian) (double *x, double*J),
		    bool (bypass_check) (double *x, double *f, double *dx),
		    int maxiter,
		    int miniter,
		    double rtol,
		    double atol,
		    double bypass_rtol,
		    double bypass_atol,
		    double residual_tol,
		    double max_dx,
		    double jmin,
		    bool random_initial,
		    performance_stat *nr_stat,
		    bool debug,
		    char *debug_file )
{
	// check load_f exist or not
	check_user_define_args( x0, load_f, load_jacobian );

	// initialize memory
	int J_size = n * n;
	int *perm = (int *) malloc ( sizeof(int) * n ); // use in lapack LU factorization
	double *x = (double *) malloc ( sizeof(double) * n );
	double *dx = (double *) malloc ( sizeof(double) * n );
	double *dx_old[2] = {0};
	double *f = (double *) malloc ( sizeof(double) * n );
	double *df = NULL;
	double *f_delta_forward = NULL;
	double *f_delta_backward = NULL;
	double *rhs = (double *) malloc ( sizeof(double) * n );
	double *J = (double *) malloc ( sizeof(double) * J_size );
	double *J_inv = NULL;
	double *J_old = NULL;
	double *D = NULL;
	double alpha;
	double beta;
	double J_norm;
	double J_inv_norm;
	double F_norm;
	double X_norm;
	bool bypass_violate;
	FILE *fout_debug = NULL;

	if ( NEWTON_DIFF_JACOBIAN != diff_type )
	{
		f_delta_forward  = (double *) malloc ( sizeof(double) * n );
		if ( diff_type == NEWTON_DIFF_CENTRAL )
		{
			f_delta_backward = (double *) malloc ( sizeof(double) * n );
		}
	}
	if ( NEWTON_BROYDEN == iterative_type )
	{
		df = (double *) malloc ( sizeof(double) * n );
		J_old = (double *) malloc ( sizeof(double) * J_size );
	}
	else if ( (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		  (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
	{
		df = (double *) malloc ( sizeof(double) * n );
	}
	if ( NEWTON_JACOBI == iterative_type )
	{
		D = (double *) malloc ( sizeof(double) * n );
	}
	if ( debug )
	{
		dx_old[0] = (double *) malloc ( sizeof(double) * n );
		dx_old[1] = (double *) malloc ( sizeof(double) * n );
		J_inv = (double *) malloc ( sizeof(double) * J_size );
	}

	newton_initialize( n, x, x0, random_initial );

	// use to output raw data for plot and analysis converge issue
	if ( debug_file )
	{
		char debug_file_name[BUFSIZ] = {0};
		sprintf( debug_file_name, "%s.raw", debug_file );
		fout_debug = fopen( debug_file_name, "w" );
		if ( !fout_debug )
		{
			fprintf( stderr, "[Error] open debug file %s fail --> %s\n", debug_file_name, strerror(errno) );
			abort();
		}
		fprintf( fout_debug, "# iter " );
		for ( int i = 0; i < n; ++i )
		{
			fprintf( fout_debug, "x%d ", i );
		}
		for ( int i = 0; i < n; ++i )
		{
			fprintf( fout_debug, "dx%d ", i );
		}
		for ( int i = 0; i < n; ++i )
		{
			fprintf( fout_debug, "f%d ", i );
		}
		fprintf( fout_debug, "\n" );
	}

	// iterative procedure
	int iter = 1;
	bool converge = false;
	double delta_ratio = 1e-9;
	double delta;
	double delta_inv;
	double x_tmp;
	while ( !converge )
	{
		if ( (iter > maxiter) && (-1 != maxiter) && (iter >= miniter) )
		{
			break;
		}

		// load RHS
		if ( (NEWTON_BROYDEN == iterative_type) ||
		     (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		     (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
		{
			memcpy( df, f, sizeof(double) * n );	
		}
		load_f( x, f );
		++(nr_stat->n_f_load);
		memcpy( rhs, f, sizeof(double) * n );	
		if ( (NEWTON_BROYDEN == iterative_type) ||
		     (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		     (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
		{
			for ( int i = 0; i < n; ++i )
			{
				df[i] = f[i] - df[i];
			}
		}
		if ( (1 == iter) && debug )
		{
			printf( "------- initial -------\n" );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.10e f=%.10e\n", i, x[i], f[i] );
			}
		}

		// construct jacobian matrix
		if ( NEWTON_CHORD_WITH_BYPASS_CHECK == iterative_type )
		{
			if ( 1 == iter )
			{
				bypass_violate = true;
			}
			else
			{
				if ( bypass_check )
				{
					bypass_violate = bypass_check( x, f, dx );
				}
				else
				{
					// user not specify bypass_check routine
					bypass_violate = __bypass_check( n, x, f, dx, bypass_rtol, bypass_atol );
				}

				if ( debug )
				{
					printf( "bypass = %s\n", bypass_violate ? "no" : "yes" );
				}
			}
		}
		if ( (NEWTON_NORMAL == iterative_type) || 
		     (NEWTON_JACOBI == iterative_type) || 
		     ((NEWTON_CHORD == iterative_type) && (1 == iter)) ||
		     ((NEWTON_CHORD_WITH_BYPASS_CHECK == iterative_type) && bypass_violate ) ||
		     ((NEWTON_BROYDEN == iterative_type) && (1 == iter)) ||
		     ((NEWTON_BROYDEN_INVERTED == iterative_type) && (1 == iter)) ||
		     ((NEWTON_BROYDEN_INVERTED_BAD == iterative_type) && (1 == iter)) 
		   )
		{
			if ( (NEWTON_DIFF_JACOBIAN  == diff_type) && load_jacobian )
			{
				// use user pre-define jacobian 
				load_jacobian( x, J );
				++(nr_stat->n_jac_load);
			}
			else
			{
				for ( int i = 0; i < n; ++i )
				{
					if ( NEWTON_DIFF_FORWARD == diff_type )
					{
						// use forward difference for better speed
						delta = x[i] * delta_ratio;
						if ( fabs(delta) < DBL_EPSILON )
						{
							if ( debug )
							{
								printf( "x[%d] is too small, delta change from %.10e to %.10e\n", i, delta, x[i] * (DBL_EPSILON / delta_ratio) );
							}
							delta = DBL_EPSILON / delta_ratio;
						}
						delta_inv = 1.0 / delta;
						x_tmp = x[i];
						x[i] += delta;
						load_f( x, f_delta_forward );
						++(nr_stat->n_f_load);
						x[i] = x_tmp;
						for ( int k = 0; k < n; ++k )
						{
							*(J + n*i + k) = (f_delta_forward[k] - f[k]) * delta_inv;
						}
					}
					else if ( NEWTON_DIFF_CENTRAL == diff_type )
					{
						// use central difference for accurate derivative approximation O(h^2)
						delta = x[i] * delta_ratio;
						if ( fabs(delta) < DBL_EPSILON )
						{
							if ( debug )
							{
								printf( "x[%d] is too small, delta change from %.10e to %.10e\n", i, delta, x[i] * (DBL_EPSILON / delta_ratio) );
							}
							delta = DBL_EPSILON / delta_ratio;
						}
						delta_inv = 1.0 / (2.0 * delta);
						x_tmp = x[i];
						x[i] += delta;
						load_f( x, f_delta_forward );
						++(nr_stat->n_f_load);
						x[i] = x_tmp;
						x[i] -= delta;
						load_f( x, f_delta_backward );
						++(nr_stat->n_f_load);
						x[i] = x_tmp;
						for ( int k = 0; k < n; ++k )
						{
							*(J + n*i + k) = (f_delta_forward[k] - f_delta_backward[k]) * delta_inv;
						}
					}
				}
			}

			// prevent matrix be singular during LU factorization, this technique will 'not' affect accuracy
			if ( jmin > 0.0 )
			{
				dense_diagonal_addition( n, J, &jmin, REAL_NUMBER );
			}

			if ( (NEWTON_BROYDEN == iterative_type) && (1 == iter) )
			{
				// save old J before factorization
				memcpy( J_old, J, sizeof(double) * J_size );	

			}
		}
		else if ( ((NEWTON_BROYDEN == iterative_type) && (iter > 1)) )
		{
			// Broyden update J
			broyden_update( n, J_old, df, dx, debug );
			memcpy( J, J_old, sizeof(double) * J_size );	
		}

		if ( debug )
		{
			if ( (NEWTON_BROYDEN_INVERTED != iterative_type) &&  
			     (NEWTON_BROYDEN_INVERTED_BAD != iterative_type) )
			{
				printf( "J = \n" );
				dense_print_matrix ( n, n, J, REAL_NUMBER );
				dense_matrix_norm ( -1, n, n, J, &J_norm, REAL_NUMBER );
				printf( "|J|_max = %.10e\n", J_norm );

				memcpy( J_inv, J, sizeof(double) * J_size );
				dense_matrix_inverse ( n, J_inv, perm, FACTOR_LU_RIGHT_LOOKING, REAL_NUMBER );
				printf( "J^-1 = \n" );
				dense_print_matrix ( n, n, J_inv, REAL_NUMBER );
				dense_matrix_norm ( -1, n, n, J_inv, &J_inv_norm, REAL_NUMBER );
				printf( "|J^-1|_max = %.10e\n", J_inv_norm );

			}
		}

		// matrix factorization A = PLU
		bool matrix_factor_ok = false;
		bool matrix_solve_ok = false;
		if ( !(NEWTON_JACOBI == iterative_type) &&
		     !((NEWTON_CHORD == iterative_type) && (iter > 1)) && 
		     !((NEWTON_CHORD_WITH_BYPASS_CHECK == iterative_type) && !bypass_violate) && 
		     !(NEWTON_BROYDEN_INVERTED == iterative_type) && 
		     !(NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
		{
			matrix_factor_ok = dense_lu_factor ( n, J, perm, FACTOR_LU_RIGHT_LOOKING, REAL_NUMBER );
			++(nr_stat->n_mat_factor);
			if ( !matrix_factor_ok )
			{
				fprintf( stderr, "[Error] LU factorization fail\n" );
				abort();
			}
		}

		// solve J(dx) = -F
		if ( debug )
		{
			// use to esitimate converge rate
			memcpy( dx_old[1], dx_old[0], sizeof(double) * n ); 
			memcpy( dx_old[0], dx, sizeof(double) * n ); 
		}
		for ( int i = 0; i < n; ++i )
		{
			rhs[i] *= -1.0;
		}
		if ( (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		     (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
		{
			if ( 1 == iter )
			{
				matrix_factor_ok = dense_matrix_inverse ( n, J, perm, FACTOR_LU_RIGHT_LOOKING, REAL_NUMBER );
				if ( !matrix_factor_ok )
				{
					fprintf( stderr, "[Error] inverse jacobian matrix fail\n" );
					abort();
				}
			}
			else 
			{
				if ( NEWTON_BROYDEN_INVERTED == iterative_type )
				{
					// Sherman-Morrison Broyden update, the same of bad Broyden but with better numerical stability
					broyden_update_sherman_morrison( n, J, df, dx, debug );
				}
				else if ( NEWTON_BROYDEN_INVERTED_BAD == iterative_type )
				{
					// bad Broyden update, more efficient but not numerical stable 
					broyden_update( n, J, dx, df, debug );
				}
			}

			if ( debug )
			{
				dense_matrix_inverse ( n, J, perm, FACTOR_LU_RIGHT_LOOKING, REAL_NUMBER );
				printf( "J = \n" );
				dense_print_matrix ( n, n, J, REAL_NUMBER );

				dense_matrix_inverse ( n, J, perm, FACTOR_LU_RIGHT_LOOKING, REAL_NUMBER );
				printf( "J^-1 = \n" );
				dense_print_matrix ( n, n, J, REAL_NUMBER );
			}

			// dx = J^-1 * -f
			alpha = 1.0;
			beta = 0.0;
			dense_matrix_vector_multiply ( n, n, &alpha, J, rhs, &beta, dx, TRANS_NONE, REAL_NUMBER );
		}
		else if ( NEWTON_JACOBI == iterative_type )
		{
			dense_matrix_get_diagonal ( n, J, D, REAL_NUMBER );
			if ( jmin > 0.0 )
			{
				for ( int i = 0; i < n; ++i )
				{
					D[i] += jmin;
				}
			}
			for ( int i = 0; i < n; ++i )
			{
				if ( 0.0 == D[i] )
				{
					fprintf( stderr, "[Error] Newton-Jacobi fail due to D[%d]=0\n", i );
					abort();
				}
				rhs[i] /= D[i];
			}
			memcpy( dx, rhs, sizeof(double) * n );
		}
		else
		{
			matrix_solve_ok = dense_solve ( n, 1, J, rhs, perm, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, REAL_NUMBER );
			++(nr_stat->n_mat_solve);
			if ( !matrix_solve_ok )
			{
				fprintf( stderr, "[Error] LU solve fail\n" );
				abort();
			}
			memcpy( dx, rhs, sizeof(double) * n );
		}

		// show solve results 
		if ( debug )
		{
			printf( "\n------- iter %d -------\n", iter );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d] new=%.10e old=%.10e dx=%.10e f=%.10e\n", i, x[i] + dx[i], x[i], dx[i], f[i] );
			}
		}

		if ( debug )
		{
			dense_vector_norm ( -1, n, dx, &X_norm, REAL_NUMBER );
			dense_vector_norm ( -1, n, f, &F_norm, REAL_NUMBER );
			printf( "[norm] |X|_max=%.10e <= |F|_max=%.10e * |J^-1|_norm=%.10e = %.10e --> %d\n", X_norm, F_norm, J_inv_norm, F_norm * J_inv_norm, F_norm * J_inv_norm > X_norm );
		}

		// modified newton 
		if ( DAMPED_DIRECT == damped_type )
		{
			// directly damped (change both direction and value)
			for ( int i = 0; i < n; ++i )
			{
				if ( fabs(dx[i]) > max_dx )
				{
					if ( debug )
					{
						printf( "[damped] change x%d from %.10e to %.10e\n", i, dx[i], ((dx[i] > 0.0) ? max_dx : -max_dx) );
					}
					dx[i] = (dx[i] > 0.0) ? max_dx : -max_dx;
				}
			}
		}

		// check stop criteria
		double diff;
		double tol;
		converge = true;
		for ( int i = 0; i < n; ++i )
		{
			// check iteration converge
			diff = fabs( dx[i] );
			tol = fabs(x[i]) * rtol + atol;
			if ( diff > tol )
			{
				if ( debug )
				{
					printf( "iter=%d x[%d]=%.10e x_new[%d]=%.10e does not converged, diff=%.10e, tol=%.10e\n", iter, i, x[i], i, x[i] + dx[i], diff, tol );
				}
				converge = false;
				break;
			}

			// check residual converge (maximum norm)
			if ( fabs(f[i]) > residual_tol )
			{
				if ( debug )
				{
					printf( "iter=%d f[%d]=%.10e does not converged residual_tol=%.10e\n", iter, i, f[i], residual_tol );
				}
				converge = false;
				break;
			}
		}

		// check converge rate
		if ( debug && (iter > 2) )
		{
			double rate;
			printf( "converge rate of each variable:\n" );
			for ( int i = 0; i < n; ++i )
			{
				rate = fabs(log( fabs(dx[i]) / fabs(dx_old[0][i]) ) / log( fabs(dx_old[0][i]) / fabs(dx_old[1][i]) ));
				printf( "x%d = %.10e (dx=%.10e dx_old0=%.10e dx_old1=%.10e\n", i, rate, dx[i], dx_old[0][i], dx_old[1][i] );
			}
		}

		// modified newton for better performance or prevent too large step
		if ( DAMPED_LINE_SEARCH == damped_type )
		{
		}

		if ( debug && fout_debug )
		{
			fprintf( fout_debug, "%d ", iter );
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "%.10e ", x[i] );
			}
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "%.10e ", dx[i] );
			}
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "%.10e ", f[i] );
			}
			fprintf( fout_debug, "\n" );
		}


		// next iteration
		for ( int i = 0; i < n; ++i )
		{
			x[i] += dx[i];
		}

		// for benchmark performance
		if ( iter < miniter )
		{
			converge = false;
		}

		++iter;
		++(nr_stat->n_iter);
	}

	// store final results
	memcpy( x_result, x, sizeof(double) * n );	
	memcpy( f_result, f, sizeof(double) * n );	

	// show newton results
	if ( debug )
	{
		printf( "\n========== Newton Converge %s in %d Iteration ==========\n", (converge ? "Success" : "Fail"), iter - 1 );
		for ( int i = 0; i < n; ++i )
		{
			printf( "x[%d]=%.10e  f=%.10e\n", i, x[i], f[i] );
		}
	}

	// release memory
	free( perm );
	free( x );
	free( dx );
	free( f );
	free( rhs );
	free( J );
	if ( diff_type != NEWTON_DIFF_JACOBIAN )
	{
		free( f_delta_forward  );
		if ( diff_type == NEWTON_DIFF_CENTRAL )
		{
			free( f_delta_backward );
		}
	}
	if ( NEWTON_BROYDEN == iterative_type )
	{
		free( df );
		free( J_old );
	}
	else if ( (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		  (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
	{
		free( df );
	}
	if ( NEWTON_JACOBI == iterative_type )
	{
		free( D );
	}
	if ( debug )
	{
		free( dx_old[0] );
		free( dx_old[1] );
	}

	return converge;
}

static void check_user_define_args ( double *x0, 
				     void (load_f) (double *x, double*f),
				     void (load_jacobian) (double *x, double*J) )
{
	if ( !x0 )
	{
		printf( "[Error] x0 is undefined\n" );
		abort();
	}
	if ( !load_f )
	{
		printf( "[Error] cannot find F(x)\n" );
		abort();
	}
	if ( !load_jacobian )
	{
		printf( "[Warning] cannot find F'(x), use approximate derivative by finite difference\n" );
	}
}

static void newton_initialize ( int n, double *x, double *x0, bool random_initial )
{
	memcpy( x, x0, sizeof(double) * n );	
}

// J = J + (df - J*dx)*dxT / |dx|^2
static void broyden_update ( int n, double *J, double *df, double *dx, bool debug )
{
	double alpha;
	double beta;
	double conjugate = false;
	double dx_square;
	double *work = (double *) malloc ( sizeof(double) * n );
	memcpy( work, df, sizeof(double) * n );

	// |dx|^2
	dense_vector_inner_product ( n, dx, dx, &dx_square, conjugate, REAL_NUMBER );

	// work = df - (J * dx)
	alpha = -1.0;
	beta = 1.0;
	dense_matrix_vector_multiply ( n, n, &alpha, J, dx, &beta, work, TRANS_NONE, REAL_NUMBER );

	// J = J + (x . yT) / |dx|^2
	alpha = 1.0 / dx_square;
	dense_maxtrix_rank_1_update ( n, n, J, &alpha, work, dx, TRANS_NONE, REAL_NUMBER );

	free( work );
}

// J = J + (dx - J*df)*dxT*J / (dxT*J*df)
static void broyden_update_sherman_morrison ( int n, double *J, double *df, double *dx, bool debug )
{
	double alpha;
	double beta;
	double conjugate = false;
	double dx_square;
	double *work1 = (double *) malloc ( sizeof(double) * n );
	double *work2 = (double *) malloc ( sizeof(double) * n );

	// work = dx - (J * df)
	memcpy( work1, dx, sizeof(double) * n );
	alpha = -1.0;
	beta = 1.0;
	dense_matrix_vector_multiply ( n, n, &alpha, J, df, &beta, work1, TRANS_NONE, REAL_NUMBER );

	// dxT*J = JT*dx
	alpha = -1.0;
	beta = 0.0;
	dense_matrix_vector_multiply ( n, n, &alpha, J, dx, &beta, work2, true, REAL_NUMBER );

	// (dxT*J) * df
	dense_vector_inner_product ( n, work2, df, &dx_square, conjugate, REAL_NUMBER );

	// J = J + (x . yT) / |dx|^2
	alpha = 1.0 / dx_square;
	dense_maxtrix_rank_1_update ( n, n, J, &alpha, work1, work2, TRANS_NONE, REAL_NUMBER );

	free( work1 );
	free( work2 );
}

static bool __bypass_check ( int n, double *x, double *f, double *dx, double bypass_rtol, double bypass_atol )
{
	bool bypass_violate = false;
	double tol;
	for ( int i = 0; i < n ; ++i )
	{
		tol = fabs(x[i] * bypass_rtol) + bypass_atol;
		if ( fabs(dx[i]) > tol )
		{
			bypass_violate = true;
			break;
		}
	}

	return bypass_violate;
}

