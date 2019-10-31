#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <errno.h>

#include "newton.h"
#include "matrix_solver.h"

// prevent double-free crash
__attribute__((always_inline)) inline void free_with_set_null ( void *ptr )
{
	if ( ptr )
	{
		free( ptr );
		ptr = NULL;
	}
}

static void check_user_define_args ( double *x0,
				     void (load_f) (double *x, double*f),
				     void (load_jacobian) (double *x, double*J) );
static void newton_initialize ( int n, double *x, double *x0, bool random_initial );

// general newton subroutines 
__attribute__((always_inline)) inline static double eval_tol( double xref, double rtol, double atol )
{
	return (fabs(xref) * rtol) + atol;
}
__attribute__((always_inline)) static double eval_local_norm( double x, double xref, double rtol, double atol )
{
	return fabs(x) / eval_tol( xref, rtol, atol );
}
static void eval_max_norm( int n, double *x, double *xref, double rtol, double atol, double *max_norm, int *max_idx );

// for broyden method, approximate jacobian
static void broyden_update ( int n, double *J, double *df, double *dx, bool debug );
static void broyden_update_sherman_morrison ( int n, double *J, double *df, double *dx, bool debug );

// for chord newton
static void chord_newton_converge_predict_iterative ( double rate, double x, double xref, double rtol, double atol, int *predict_iter, double *predict_norm );
static void chord_newton_converge_predict_approximate ( double rate, double norm, int *predict_iter, double *predict_norm );
static bool __bypass_check ( int n, double *x, double *f, double *dx, double bypass_rtol, double bypass_atol );

bool newton_solve ( newton_iterative_type iterative_type, 
		    newton_damped_type damped_type,
		    newton_rescue_type rescue_type,
		    newton_derivative_type diff_type,
		    int n,
		    double *x0,
		    double *x_ans,
		    double *x_result,
		    double *f_result,
		    void (load_f) (double *x, double*f),
		    void (load_jacobian) (double *x, double*J),
		    bool (bypass_check) (double *x, double *f, double *dx),
		    int maxiter,
		    int miniter,
		    double delta_rtol,
		    double delta_atol,
		    double residual_rtol,
		    double residual_atol,
		    double bypass_rtol,
		    double bypass_atol,
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
	double *dx2 = NULL;
	double *f = (double *) malloc ( sizeof(double) * n );
	double *f_old = (double *) malloc ( sizeof(double) * n );
	double *df = (double *) malloc ( sizeof(double) * n );
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

	// ---------------------------------
	// necessarily memory allocation with different method
	// ---------------------------------
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
		J_old = (double *) malloc ( sizeof(double) * J_size );
	}
	else if ( (NEWTON_CHORD == iterative_type) ||
		  (NEWTON_CHORD_WITH_BYPASS_CHECK == iterative_type) )
	{
		dx_old[0] = (double *) malloc ( sizeof(double) * n );
		dx_old[1] = (double *) malloc ( sizeof(double) * n );
	}
	else if ( NEWTON_JACOBI == iterative_type )
	{
		D = (double *) malloc ( sizeof(double) * n );
	}

	// ---------------------------------
	// debug used memory allocation
	// ---------------------------------
	if ( debug )
	{
		dx_old[0] = (double *) malloc ( sizeof(double) * n );
		dx_old[1] = (double *) malloc ( sizeof(double) * n );
		dx2 = (double *) malloc ( sizeof(double) * n );
		J_inv = (double *) malloc ( sizeof(double) * J_size );
	}

	newton_initialize( n, x, x0, random_initial );

	// ---------------------------------
	// use to output raw data for plot and analysis converge issue
	// ---------------------------------
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

	// ---------------------------------
	// iterative procedure
	// ---------------------------------
	int iter = 1;
	int idx;
	int max_dx_idx = -1;
	int max_f_idx = -1;
	double dx_max_norm = -1;
	double f_max_norm = -1;
	double local_norm;
	double delta_ratio = 1e-9;
	double delta;
	double delta_inv;
	double diff;
	double tol;
	double x_tmp;
	bool delta_converge = false;
	bool f_converge = false;
	bool nr_converge = false;
	while ( !nr_converge )
	{
		if ( (iter > maxiter) && (-1 != maxiter) && (iter >= miniter) )
		{
			break;
		}

		// ---------------------------------
		// load RHS
		// ---------------------------------
		memcpy( f_old, f, sizeof(double) * n );	
		load_f( x, f );
		++(nr_stat->n_f_load);
		memcpy( rhs, f, sizeof(double) * n );	
		for ( int i = 0; i < n; ++i )
		{
			df[i] = f[i] - f_old[i];
		}
		if ( debug )
		{
			printf( "\n------- iter %d -------\n", iter );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.15le f[%d]=%.15le\n", i, x[i], i, f[i] );
			}
		}


		// check residue converged after load f
		// ‖.‖≡ ‖.‖∞
		// ‖f‖ = ‖f/tol‖∞
		if ( !f_converge )
		{
			// check delta converge
			// ‖.‖≡ ‖.‖∞
			// ‖f‖ = ‖f/tol‖∞
			eval_max_norm( n, f, f, residual_rtol, residual_atol, &f_max_norm, &max_f_idx );
			if ( f_max_norm > 1 )
			{
				f_converge = false;
				if ( debug )
				{
					idx = max_f_idx;
					tol = eval_tol( f[idx], residual_rtol, residual_atol );
					printf( "iter=%d norm=%.15le f[%d]=%.15le f_old[%d]=%.15le residue non-converged, df=%.15le, tol=%.15le\n", iter, f_max_norm, idx, f[idx], idx, f[idx] - df[idx], df[idx], tol );
				}
			}
			else
			{
				f_converge = true;
			}
		}

		// complete newton if both residue and delta converge 
		nr_converge = (delta_converge && f_converge);
		if ( nr_converge )
		{
			if ( debug )
			{
				printf( "[converge info] iter=%d both delta and f converge, skip load jacobian and matrix solve\n", iter );
			}
			break;
		}

		// ---------------------------------
		// construct jacobian matrix
		// ---------------------------------
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
								printf( "x[%d] is too small, delta change from %.15le to %.15le\n", i, delta, x[i] * (DBL_EPSILON / delta_ratio) );
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
								printf( "x[%d] is too small, delta change from %.15le to %.15le\n", i, delta, x[i] * (DBL_EPSILON / delta_ratio) );
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
				printf( "|J|_max = %.15le\n", J_norm );

				memcpy( J_inv, J, sizeof(double) * J_size );
				dense_matrix_inverse ( n, J_inv, perm, FACTOR_LU_RIGHT_LOOKING, REAL_NUMBER );
				printf( "J^-1 = \n" );
				dense_print_matrix ( n, n, J_inv, REAL_NUMBER );
				dense_matrix_norm ( -1, n, n, J_inv, &J_inv_norm, REAL_NUMBER );
				printf( "|J^-1|_max = %.15le\n", J_inv_norm );

			}
		}

		// ---------------------------------
		// matrix factorization A = PLU
		// ---------------------------------
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

		// ---------------------------------
		// matix solve J*dx = -F
		// ---------------------------------
		if ( debug )
		{
			// use to estimate converge rate
			memcpy( dx_old[1], dx_old[0], sizeof(double) * n ); 
			memcpy( dx_old[0], dx, sizeof(double) * n ); 

			if ( 2 == iter ) // use for chord newton rate estimate
			{
				memcpy( dx2, dx, sizeof(double) * n ); 
			}
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

			// dx = J⁻¹ * -f
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
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d] new=%.15le old=%.15le dx=%.15le f=%.15le\n", i, x[i] + dx[i], x[i], dx[i], f[i] );
			}
		}

		if ( debug )
		{
			dense_vector_norm ( -1, n, dx, &X_norm, REAL_NUMBER );
			dense_vector_norm ( -1, n, f, &F_norm, REAL_NUMBER );
			printf( "[norm] |X|_max=%.15le <= |F|_max=%.15le * |J^-1|_norm=%.15le = %.15le --> %d\n", X_norm, F_norm, J_inv_norm, F_norm * J_inv_norm, F_norm * J_inv_norm > X_norm );
		}

		// ---------------------------------
		// damped newton, limit dx
		// ---------------------------------
		if ( DAMPED_DIRECT == damped_type )
		{
			// directly damped (change both direction and value)
			for ( int i = 0; i < n; ++i )
			{
				if ( fabs(dx[i]) > max_dx )
				{
					if ( debug )
					{
						printf( "[damped] change x%d from %.15le to %.15le\n", i, dx[i], ((dx[i] > 0.0) ? max_dx : -max_dx) );
					}
					dx[i] = (dx[i] > 0.0) ? max_dx : -max_dx;
				}
			}
		}

		// ---------------------------------
		// check delta converge after matrix solve
		// ---------------------------------
		if ( !delta_converge )
		{
			// check delta converge
			// ‖.‖≡ ‖.‖∞
			// ‖dx‖ = ‖dx/tol‖∞
			eval_max_norm( n, dx, x, delta_rtol, delta_atol, &dx_max_norm, &max_dx_idx );
			if ( dx_max_norm > 1 )
			{
				delta_converge = false;
				if ( debug )
				{
					idx = max_dx_idx;
					tol = eval_tol( x[idx], residual_rtol, residual_atol );
					printf( "iter=%d norm=%.15le x[%d]=%.15le x_new[%d]=%.15le delta non-converged, dx=%.15le, tol=%.15le\n", iter, dx_max_norm, idx, x[idx], idx, x[idx] + dx[idx], dx[idx], tol );
				}
			}
			else
			{
				delta_converge = true;
			}
		}

		nr_converge = (delta_converge && f_converge);
		if ( debug )
		{
			printf( "iter=%d nr_converge=%d delta_converge=%d f_converge=%d\n", iter, nr_converge, delta_converge, f_converge );
		}

		// ---------------------------------
		// check converge order and chord newton linear rate
		// ---------------------------------
		if ( debug && (iter > 2) )
		{
			double converge_order;
			printf( "converge order of each variable:\n" );
			for ( int i = 0; i < n; ++i )
			{
				converge_order = fabs(log( fabs(dx[i]) / fabs(dx_old[0][i]) ) / log( fabs(dx_old[0][i]) / fabs(dx_old[1][i]) ));
				printf( "x%d = %.15le (dx=%.15le dx_old0=%.15le dx_old1=%.15le\n", i, converge_order, dx[i], dx_old[0][i], dx_old[1][i] );
			}
		}

		if ( debug && (iter > 1) )
		{
			if ( NEWTON_CHORD == iterative_type )
			{
				// estimate linear rate
				double rate_dx = NAN;
				double rate_dx_avg;
				double fixed_point_dist_estimate;
				double fixed_point_dist_estimate_avg;
				double fixed_point_dist_exact;
				printf( "linear converge rate_dx and fixed-point distance of chord newton:\n" );
				for ( int i = 0; i < n; ++i )
				{
					// ‖eₖ‖ = ‖x* - xₖ‖
					// ‖eₖ₊₁‖ ≈ |ρ| * ‖eₖ‖
					rate_dx = fabs(dx[i]) / fabs(dx_old[0][i]);
					rate_dx_avg = exp( (1.0/(iter-2)) * log(fabs(dx[i]) / fabs(dx2[i])) );
					fixed_point_dist_estimate = fabs((rate_dx / (1 - rate_dx)) * dx[i]);
					fixed_point_dist_estimate_avg = fabs((rate_dx_avg / (1 - rate_dx_avg)) * dx[i]);
					if ( NULL == x_ans )
					{
						printf( "rate_dx%d = %.15le, rate_dx_avg = %.15le, dist_estimate = %.15le, dist_estimate_avg = %.15le\n", i, rate_dx, rate_dx_avg, fixed_point_dist_estimate, fixed_point_dist_estimate_avg );
					}
					else
					{ fixed_point_dist_exact = fabs((x[i] + dx[i]) - x_ans[i]);
						printf( "rate_dx%d = %.15le, rate_dx_avg = %.15le, dist_estimate= %.15le, dist_estimate_avg = %.15le, dist_exact = %.15le\n", i, rate_dx, rate_dx_avg, fixed_point_dist_estimate, fixed_point_dist_estimate_avg, fixed_point_dist_exact );
					}
				}

				double rate_f = NAN;
				for ( int i = 0; i < n; ++i )
				{
					rate_f = fabs(f[i]) / fabs(f_old[i]);
					printf( "rate_f%d = %.15le\n", i, rate_f );
				}

				// estimate need how many following iterations for statisfy converge criteria
				int n_iter_for_delta_converge;
				double dx_predict_norm;
				int n_iter_for_f_converge;
				double f_predict_norm;
				if ( !delta_converge )
				{
					// ‖dx‖ * |rate|ⁿ ≤ 1
					chord_newton_converge_predict_iterative( rate_dx, dx[max_dx_idx], x[max_dx_idx], delta_rtol, delta_atol, &n_iter_for_delta_converge, &dx_predict_norm );
					printf( "[converge predict iterative]   need %d iter for delta norm converge to %.15le\n", n_iter_for_delta_converge, dx_predict_norm );
					chord_newton_converge_predict_approximate( rate_dx, dx_max_norm, &n_iter_for_delta_converge, &dx_predict_norm );
					printf( "[converge predict approximate] need %d iter for delta norm converge to %.15le\n", n_iter_for_delta_converge, dx_predict_norm );
				}
				if ( !f_converge )
				{
					// ‖f‖ * |rate|ⁿ ≤ 1
					chord_newton_converge_predict_iterative( rate_f, f[max_f_idx], NAN, residual_rtol, residual_atol, &n_iter_for_f_converge, &f_predict_norm );
					printf( "[converge predict iterative]   need %d iter for f norm converge to %.15le\n", n_iter_for_f_converge, f_predict_norm );
					chord_newton_converge_predict_approximate( rate_f, f_max_norm, &n_iter_for_f_converge, &f_predict_norm );
					printf( "[converge predict approximate] need %d iter for f norm converge to %.15le\n", n_iter_for_f_converge, f_predict_norm );
				}
			}
		}

		// ---------------------------------
		// modified newton for better performance or prevent too large step
		// ---------------------------------
		if ( DAMPED_LINE_SEARCH == damped_type )
		{
		}

		if ( debug && fout_debug )
		{
			fprintf( fout_debug, "%d ", iter );
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "%.15le ", x[i] );
			}
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "%.15le ", dx[i] );
			}
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "%.15le ", f[i] );
			}
			fprintf( fout_debug, "\n" );
		}


		// ---------------------------------
		// next iteration
		// ---------------------------------
		for ( int i = 0; i < n; ++i )
		{
			x[i] += dx[i];
		}

		// for benchmark performance
		if ( iter < miniter )
		{
			nr_converge = false;
		}

		if ( !nr_converge )
		{
			printf( "[converge info] iter=%-3d dx_norm=%.15le (id=%-3d), f_norm=%.15le (id=%-3d)\n", iter, dx_max_norm, max_dx_idx, f_max_norm, max_f_idx );
		}

		++iter;
		++(nr_stat->n_iter);
	}

	// ---------------------------------
	// complete newton iterations, store final results
	// ---------------------------------
	memcpy( x_result, x, sizeof(double) * n );	
	memcpy( f_result, f, sizeof(double) * n );	

	// show newton results
	if ( debug )
	{
		printf( "\n========== Newton Converge %s in %d Iteration ==========\n", (nr_converge ? "Success" : "Fail"), iter );
		for ( int i = 0; i < n; ++i )
		{
			printf( "x[%d]=%.15le  f=%.15le\n", i, x[i], f[i] );
		}
		printf( "[norm] dx_norm=%.15le (id=%d), f_norm=%.15le (id=%d)\n", dx_max_norm, max_dx_idx, f_max_norm, max_f_idx );
	}

	// ---------------------------------
	// release memory
	// ---------------------------------
	free_with_set_null( perm );
	free_with_set_null( x );
	free_with_set_null( dx );
	free_with_set_null( f );
	free_with_set_null( f_old );
	free_with_set_null( df );
	free_with_set_null( rhs );
	free_with_set_null( J );
	if ( diff_type != NEWTON_DIFF_JACOBIAN )
	{
		free_with_set_null( f_delta_forward  );
		if ( diff_type == NEWTON_DIFF_CENTRAL )
		{
			free_with_set_null( f_delta_backward );
		}
	}

	if ( NEWTON_BROYDEN == iterative_type )
	{
		free_with_set_null( J_old );
	}
	else if ( (NEWTON_CHORD == iterative_type) ||
		  (NEWTON_CHORD_WITH_BYPASS_CHECK == iterative_type) )
	{
		free_with_set_null( dx_old[0] );
		free_with_set_null( dx_old[1] );
	}
	else if ( NEWTON_JACOBI == iterative_type )
	{
		free_with_set_null( D );
	}

	if ( debug )
	{
		free_with_set_null( dx_old[0] );
		free_with_set_null( dx_old[1] );
		free_with_set_null( dx2 );
	}

	return nr_converge;
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

static void eval_max_norm( int n, double *x, double *xref, double rtol, double atol, double *max_norm, int *max_idx )
{
	double local_norm;
	double _max_norm = -1;
	int _max_idx = -1;
	for ( int i = 0; i < n; ++i )
	{
		local_norm = eval_local_norm( x[i], xref[i], rtol, atol );
		if ( local_norm > _max_norm )
		{
			_max_norm = local_norm;
			_max_idx = i;
		}
	}
	*max_norm = _max_norm;
	*max_idx = _max_idx;
}

// J = J + ((Δf - J*Δx)*Δxᵀ) / ‖Δx‖²
static void broyden_update ( int n, double *J, double *df, double *dx, bool debug )
{
	double alpha;
	double beta;
	double conjugate = false;
	double dx_square;
	double *work = (double *) malloc ( sizeof(double) * n );
	memcpy( work, df, sizeof(double) * n );

	// ‖Δx‖²
	dense_vector_inner_product ( n, dx, dx, &dx_square, conjugate, REAL_NUMBER );

	// work = (Δf - J*Δx)
	alpha = -1.0;
	beta = 1.0;
	dense_matrix_vector_multiply ( n, n, &alpha, J, dx, &beta, work, TRANS_NONE, REAL_NUMBER );

	// J = J + (work * Δxᵀ) / ‖Δx‖²
	alpha = 1.0 / dx_square;
	dense_maxtrix_rank_1_update ( n, n, J, &alpha, work, dx, TRANS_NONE, REAL_NUMBER );

	free( work );
}

// better numerical stability then bad broyden inverted method
// J⁻¹ = J⁻¹ + ((Δx - J⁻¹*Δf)*Δxᵀ*J⁻¹) / (Δxᵀ*J⁻¹*Δf)
static void broyden_update_sherman_morrison ( int n, double *J, double *df, double *dx, bool debug )
{
	double alpha;
	double beta;
	double conjugate = false;
	double dx_square;
	double *work1 = (double *) malloc ( sizeof(double) * n );
	double *work2 = (double *) malloc ( sizeof(double) * n );

	// work1 = (Δx - J⁻¹*Δf)
	memcpy( work1, dx, sizeof(double) * n );
	alpha = -1.0;
	beta = 1.0;
	dense_matrix_vector_multiply ( n, n, &alpha, J, df, &beta, work1, TRANS_NONE, REAL_NUMBER );

	// work2 = Δxᵀ*J⁻¹ = (J⁻¹ᵀ)*Δx
	alpha = -1.0;
	beta = 0.0;
	dense_matrix_vector_multiply ( n, n, &alpha, J, dx, &beta, work2, true, REAL_NUMBER );

	// ‖Δx‖² = (Δxᵀ*J⁻¹*Δf) = work2 * Δf
	dense_vector_inner_product ( n, work2, df, &dx_square, conjugate, REAL_NUMBER );

	// J⁻¹ = J⁻¹ + (work1⨯ work2) / ‖Δx‖² 
	alpha = 1.0 / dx_square;
	dense_maxtrix_rank_1_update ( n, n, J, &alpha, work1, work2, TRANS_NONE, REAL_NUMBER );

	free( work1 );
	free( work2 );
}

static void chord_newton_converge_predict_iterative ( double rate, double x, double xref, double rtol, double atol, int *predict_iter, double *predict_norm )
{
	double norm;
	int iter = 0;
	bool xref_use_x = false;

	if ( rate >= 1 )
	{
		*predict_iter = 0xffffffff; // INT_MAX
		*predict_norm = NAN; 
		return;
	}

	if ( isnan( xref ) )
	{
		xref_use_x = true;	
	}

	if ( xref_use_x )
	{
		xref = x;
	}
	norm = eval_local_norm( x, xref, rtol, atol );
	while ( norm > 1 )
	{
		x *= rate;
		++iter;
		if ( xref_use_x )
		{
			xref = x;
		}
		norm = eval_local_norm( x, xref, rtol, atol );
	}

	*predict_iter = iter;
	*predict_norm = norm;
}

static void chord_newton_converge_predict_approximate ( double rate, double norm, int *predict_iter, double *predict_norm )
{
	// ‖f‖ * |rate|ⁿ ≤ 1
	*predict_iter = (int) ceil( log(1.0 / norm) / log(rate) );
	*predict_norm = norm * exp( *predict_iter * log( rate ) );
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
