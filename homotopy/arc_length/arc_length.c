#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <errno.h>

#include "homotopy.h"
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

// load arc length constrain
static double eval_arc_length_constrain ( int n, double *x, double *xc, double p, double pc, double dt );

// general newton subroutines 
__attribute__((always_inline)) inline static double eval_tol( double xref, double rtol, double atol )
{
	return (fabs(xref) * rtol) + atol;
}
__attribute__((always_inline)) inline static double eval_local_norm( double x, double xref, double rtol, double atol )
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

// for line-search modify newton
static double eval_f_optimization ( void (load_f) (double *x, double*f), double *x, double *dx, double a, newton_param_t *newton_param );
static double interval_halving ( void (load_f) (double *v, double*f), double *v, double *dv, double a, double b, newton_param_t *newton_param );
static double golden_section ( void (load_f) (double *v, double*f), double *v, double *dv, double a, double b, newton_param_t *newton_param );

// dump debug
static void dump_debug_data ( FILE *fout, int n, int iter, double *x, double *dx, double *f );
bool arc_length_bbd_newton_solve ( 
		newton_param_t *newton_param,
		int *perm, // permuation for matrix ordering
		double *J, // jacobian
		double *pp,
		double *x0, // initial x
	        double *x_ans, // user give solution x*
		double p0, // initial p
		double *x_result, // final x
		double *p_result, // final x
		double *f_result, // final f(x)
		double *g_result, // final g(x)
		double *xc, // hyper circle central (xc,pc)
		double pc, // hyper circle central (xc,pc)
		double dt, // arc length
		void (load_f) (double *x, double*f),
		void (load_df_dp) (double *x, double*f),
		void (load_jacobian) (double *x, double*J),
		char *debug_file )
{
	// assign newton parameters to local
	newton_iterative_type iterative_type = newton_param->iterative_type;
	newton_damped_type damped_type = newton_param->damped_type;
	newton_rescue_type rescue_type = newton_param->rescue_type;
	newton_derivative_type diff_type = newton_param->diff_type;
	int n = newton_param->n;
	int maxiter = newton_param->maxiter;
	int miniter = newton_param->miniter;
	double delta_rtol = newton_param->delta_rtol;
	double delta_atol = newton_param->delta_atol;
	double residual_rtol = newton_param->residual_rtol;
	double residual_atol = newton_param->residual_atol;
	double max_dx = newton_param->max_dx;
	double jmin = newton_param->jmin;
	double dp;
	double p;
	double g;
	bool debug = newton_param->debug;
	bool random_initial = newton_param->random_initial;
	performance_stat_t *nr_stat = &(newton_param->nr_stat);

	// check load_f exist or not
	check_user_define_args( x0, load_f, load_jacobian );

	// initialize memory
	int J_size = n * n;
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
	double *J_inv = NULL;
	double *J_old = NULL;
	double *D = NULL;
	double alpha;
	double beta;
	double J_norm;
	double J_inv_norm;
	double F_norm;
	double X_norm;
	FILE *fout_debug = NULL;

	// BBD 
	double *bbd_tmp_rhs = (double *) malloc ( sizeof(double) * n );
	double bbd_ext_J;
	double bbd_ext_f;
	double *A21 = (double *) malloc ( sizeof(double) * n );
	double *A12 = (double *) malloc ( sizeof(double) * n );
	double A22;

	if ( NULL == J )
	{
		J = (double *) malloc ( sizeof(double) * J_size );
	}
	if ( NULL == perm )
	{
		perm = (int *) malloc ( sizeof(int) * n ); // use in lapack LU factorization
	}

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
	p = p0;
	*pp = p;
	if ( debug )
	{
		if ( dt > 0 )
		{
			printf( "dt=%.15le\n", dt );
			for ( int i = 0; i < n; ++i )
			{
				printf( "xc[%d]=%.15le\n", i, xc[i] );
			}
			printf( "pc=%.15le", pc );
		}
	}

	// ---------------------------------
	// use to output raw data for plot and analysis converge issue
	// ---------------------------------
	if ( debug && debug_file )
	{
		char debug_file_name[BUFSIZ] = {0};
		sprintf( debug_file_name, "%s.nr", debug_file );
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
	int iter = 0;
	int idx;
	int max_dx_idx = -1;
	int max_f_idx = -1;
	double dx_max_norm = NAN;
	double f_max_norm = NAN;
	double dp_norm = NAN;
	double g_norm = NAN;
	double local_norm;
	double delta_ratio = 1e-9;
	double delta;
	double delta_inv;
	double diff;
	double tol;
	double x_tmp;
	bool delta_converge = false;
	bool p_converge = false;
	bool f_converge = false;
	bool g_converge = false;
	bool nr_converge = false;
	while ( true )
	{
		++iter;
		++(nr_stat->n_iter);

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

		// load arc length constrain
		if ( dt > 0 )
		{
			g = eval_arc_length_constrain( n, x, xc, p, pc, dt );
		}

		if ( debug )
		{
			printf( "\n------- iter %d -------\n", iter );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.15le f[%d]=%.15le\n", i, x[i], i, f[i] );
			}
			if ( dt > 0 )
			{
				printf( "p=%.15le g=%.15le\n", p, g );
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
		if ( !g_converge )
		{
			g_norm = eval_local_norm( g, g, residual_rtol, residual_atol );
			if ( g_norm > 1 )
			{
				g_converge = false;
				if ( debug )
				{
					printf( "iter=%d g non-converged, norm=%.15le\n", iter, g_norm );
				}
			}
			else
			{
				g_converge = true;
			}
		}

		// complete newton if both residue and delta converge 
		nr_converge = (delta_converge && f_converge && p_converge && g_converge);
		if ( nr_converge && (iter >= miniter) )
		{
			if ( debug )
			{
				printf( "[converge info] iter=%d both delta and f converge, skip load jacobian and matrix solve\n", iter );
			}
			dump_debug_data( fout_debug, n, iter, x, dx, f );
			break;
		}

		// ---------------------------------
		// construct jacobian matrix
		// ---------------------------------
		if ( (NEWTON_NORMAL == iterative_type) || 
		     (NEWTON_JACOBI == iterative_type) || 
		     ((NEWTON_CHORD == iterative_type) && (1 == iter)) ||
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
		if ( dt > 0 )
		{
			// ---------------------------------
			// BBD solve
			// A₁₁ = ∂f/∂x = J(x)
			// A₁₂ = ∂f/∂p
			// A₂₁ = ∂g/∂x
			// A₂₂ = ∂g/∂p
			// Δp = 1/(A₂₂ - A₂₁*A₁₁⁻¹*A₁₂)*(-g + A₂₁*A₁₁⁻¹*f)
			// Δx = A₁₁⁻¹*(-f - A₁₂*Δp)
			// g(x,p) = (‖x-xc‖² + ‖p-pc‖²) - Δt² 
			// ---------------------------------

			// solve Δp = 1/(A₂₂  - A₂₁*A₁₁⁻¹*A₁₂)*(-g + A₂₁*A₁₁⁻¹*f)
			load_df_dp ( x, A12 ); 
			for ( int i = 0; i < n; ++i )
			{
				A21[i] = 2 * (x[i] - xc[i]);
			}
			A22 = 2 * (p - pc);
			memcpy( bbd_tmp_rhs, A12, sizeof(double) * n ); 
			matrix_solve_ok = dense_solve ( n, 1, J, bbd_tmp_rhs, perm, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, REAL_NUMBER );
			++(nr_stat->n_mat_solve);
			if ( !matrix_solve_ok )
			{
				fprintf( stderr, "[Error] LU solve fail\n" );
				abort();
			}
			bbd_ext_J = 0;
			for ( int i = 0; i < n; ++i )
			{
				bbd_ext_J += A21[i] * bbd_tmp_rhs[i];
			}
			bbd_ext_J = A22 - bbd_ext_J;
			memcpy( bbd_tmp_rhs, f, sizeof(double) * n ); 
			matrix_solve_ok = dense_solve ( n, 1, J, bbd_tmp_rhs, perm, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, REAL_NUMBER );
			++(nr_stat->n_mat_solve);
			if ( !matrix_solve_ok )
			{
				fprintf( stderr, "[Error] LU solve fail\n" );
				abort();
			}
			bbd_ext_f = 0;
			for ( int i = 0; i < n; ++i )
			{
				bbd_ext_f += A21[i] * bbd_tmp_rhs[i];
			}
			bbd_ext_f = -g + bbd_ext_f;
			dp = bbd_ext_f / bbd_ext_J;
			if ( debug )
			{
				printf( "==================== BBD ====================\n" );
				printf( "J_external=%.15le f_external=%.15le g=%.15le\n", bbd_ext_J, bbd_ext_f, g );
				printf( "Δp=%.15le pₖ₊₁=%.15le pₖ=%.15le\n", dp, p + dp, p );
				printf( "A12=" );
				for ( int i = 0; i < n; ++i )
				{
					printf( " %.15le", A12[i] );
				}
				printf( "\n" );
				printf( "A21=" );
				for ( int i = 0; i < n; ++i )
				{
					printf( " %.15le", A21[i] );
				}
				printf( "\n" );
				printf( "A22= %.15le\n", A22 );
				printf( "=============================================\n" );
			}
			p += dp;
			*pp = p;

			// solve Δx = A₁₁⁻¹*(-f - A₁₂*Δp)
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
				bbd_tmp_rhs[i] = -f[i] - A12[i] * dp;
			}
			matrix_solve_ok = dense_solve ( n, 1, J, bbd_tmp_rhs, perm, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, REAL_NUMBER );
			++(nr_stat->n_mat_solve);
			if ( !matrix_solve_ok )
			{
				fprintf( stderr, "[Error] LU solve fail\n" );
				abort();
			}
			memcpy( dx, bbd_tmp_rhs, sizeof(double) * n );
			if ( debug )
			{
				for ( int i = 0; i < n; ++i )
				{
					printf( "x[%d] new=%.15le old=%.15le dx=%.15le f=%.15le\n", i, x[i] + dx[i], x[i], dx[i], f[i] );
				}
				dense_vector_norm ( -1, n, dx, &X_norm, REAL_NUMBER );
				dense_vector_norm ( -1, n, f, &F_norm, REAL_NUMBER );
				printf( "[norm] |X|_max=%.15le <= |F|_max=%.15le * |J^-1|_norm=%.15le = %.15le --> %d\n", X_norm, F_norm, J_inv_norm, F_norm * J_inv_norm, F_norm * J_inv_norm > X_norm );
			}
		}
		else
		{
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
		if ( !p_converge)
		{
			dp_norm = eval_local_norm( dp, p, delta_rtol, delta_atol );
			if ( dp_norm > 1 )
			{
				p_converge = false;
				if ( debug )
				{
					printf( "iter=%d dp non-converged, norm=%.15le\n", iter, dp_norm );
				}
			}
			else
			{
				p_converge = true;
			}
		}

		nr_converge = (delta_converge && f_converge && p_converge && g_converge);

		if ( nr_converge && (iter >= miniter) )
		{
			printf( "[converge info] iter=%d both delta and f converge\n", iter );
			dump_debug_data( fout_debug, n, iter, x, dx, f );
			break;
		}
		else
		{
			printf( "[converge info] iter=%-3d dx_norm=%.15le (id=%-3d), f_norm=%.15le (id=%-3d)\n", iter, dx_max_norm, max_dx_idx, f_max_norm, max_f_idx );
			if ( debug )
			{
				printf( "iter=%d nr_converge=%d delta_converge=%d f_converge=%d\n", iter, nr_converge, delta_converge, f_converge );
			}
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
					{ 
						fixed_point_dist_exact = fabs((x[i] + dx[i]) - x_ans[i]);
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
			// xₖ₊₁ = xₖ + 𝛼*dx, line search 𝛼 cause minimum ‖f(xₖ₊₁)‖ 
			//double optima_a = interval_halving ( load_f, x, dx, 0, 1, newton_param );
			double optima_a = golden_section ( load_f, x, dx, 0, 1, newton_param );
			if ( optima_a < 1 )
			{
				if ( newton_param->debug )
				{
					printf( "[line_search] optimized a=%.10le\n", optima_a );
				}

				// scale dx
				for ( int i = 0; i < n; ++i )
				{
					if ( newton_param->debug )
					{
						printf( "[damped newton] modify dx[%d] %.10le -> %.10le\n", i, dx[i], optima_a * dx[i]);
					}
					dx[i] *= optima_a;
				}
			}
		}
		

		// ---------------------------------
		// next iteration
		// ---------------------------------
		dump_debug_data( fout_debug, n, iter, x, dx, f );
		for ( int i = 0; i < n; ++i )
		{
			x[i] += dx[i];
		}
	}

	printf( "[converge info] iter=%-3d dx_norm=%.15le (id=%-3d), f_norm=%.15le (id=%-3d)\n", iter, dx_max_norm, max_dx_idx, f_max_norm, max_f_idx );

	// ---------------------------------
	// complete newton iterations, store final results
	// ---------------------------------
	memcpy( x_result, x, sizeof(double) * n );	
	memcpy( f_result, f, sizeof(double) * n );	
	if ( dt > 0 )
	{
		*p_result = p;
		*g_result = g;
	}

	// show newton results
	if ( debug )
	{
		printf( "\n========== Newton Converge %s in %d Iteration ==========\n", (nr_converge ? "Success" : "Fail"), iter );
		for ( int i = 0; i < n; ++i )
		{
			printf( "x[%d]=%.15le  f=%.15le\n", i, x[i], f[i] );
		}
		printf( "[norm] dx_norm=%.15le (id=%d), f_norm=%.15le (id=%d)\n", dx_max_norm, max_dx_idx, f_max_norm, max_f_idx );
		printf( "[norm] dp_norm=%.15le, g_norm=%.15le\n", dp_norm, g_norm );
	}

	// ---------------------------------
	// release memory
	// ---------------------------------
	free_with_set_null( x );
	free_with_set_null( dx );
	free_with_set_null( f );
	free_with_set_null( f_old );
	free_with_set_null( df );
	free_with_set_null( rhs );
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

static double eval_f_optimization ( void (load_f) (double *x, double*f), double *x, double *dx, double a, newton_param_t *newton_param )
{
	static double *x_buf = NULL;
	static double *f_buf = NULL;
	int n = newton_param->n;
	if ( !x_buf )
	{
		x_buf = malloc( sizeof(double) * n );
		f_buf = malloc( sizeof(double) * n );
	}

	double f_max_norm = -DBL_MAX;
	int max_f_idx = -1;
	for ( int i = 0; i < n; ++i )
	{
		x_buf[i] = x[i] + (a * dx[i]);
	}
	load_f( x_buf, f_buf );
	++(newton_param->nr_stat.n_f_load);
	//eval_max_norm( n, f_buf, f_buf, newton_param->residual_rtol, newton_param->residual_atol, &f_max_norm, &max_f_idx );
	for ( int i = 0; i < n; ++i )
	{
		if ( fabs(f_buf[i]) > f_max_norm )
		{
			f_max_norm = fabs(f_buf[i]);
			max_f_idx = i;
		}
		
	}
	return f_max_norm;
}

// <------  L  ----->
// a___x₁______x₂___b
// 
// assume problem is convex hull, then
// x₁ = a + 0.25*L
// x₂ = b - 0.25*L
// if f(x₁) < f(x₂) then prun x₂ ~ b
// if f(x₂) < f(x₁) then prun a ~ x₁ 
// if f(x₁) = f(x₂) then prun a ~ x₁ and x₂ ~ b
// 
// each iteration call two f(x)
// linear converge rate 0.75
//
static double interval_halving ( void (load_f) (double *v, double*f), double *v, double *dv, double a, double b, newton_param_t *newton_param )
{
	int iter;
	double x_optima = NAN;
	double f_optima = NAN;
	double f1;
	double f2;
	double x1;
	double x2;
	double x_origin;
	double f_origin;
	double l;
	double delta;
	double halving_ratio = 0.25;
	double tol = newton_param->line_search_tol;

	x_origin = b;
	f_origin = eval_f_optimization ( load_f, v, dv, b, newton_param );

	iter = 1;
	l = b - a;
	delta = l * halving_ratio;
	x1 = a + delta;
	x2 = b - delta;
	f1 = eval_f_optimization ( load_f, v, dv, x1, newton_param );
	f2 = eval_f_optimization ( load_f, v, dv, x2, newton_param );

	if ( f1 <= f2 )
	{
		x_optima = x1;
		f_optima = f1;
	}
	else if ( f2 < f1 )
	{
		x_optima = x2;
		f_optima = f2;
	}

	while( l > tol ) 
	{
		if ( newton_param->debug )
		{
			printf( "[line_search] i=%d l=%.10le a=%.10le b=%.10le x1=%.10le x2=%.10le f1=%.10le f2=%.10le f_optima=%.10le\n", iter, l, a, b, x1, x2, f1, f2, f_optima );
		}

		if ( f1 < f2 )
		{
			if ( f1 < f_optima )
			{
				x_optima = x1;
				f_optima = f1;
			}
			b = x2;
		}
		else if ( f2 < f1 )
		{
			if ( f2 < f_optima )
			{
				x_optima = x2;
				f_optima = f2;
			}
			a = x1;
		}
		else
		{
			a = x1;
			b = x2;
		}

		++iter;
		l = b - a;
		delta = l * halving_ratio;
		x1 = a + delta;
		x2 = b - delta;
		f1 = eval_f_optimization ( load_f, v, dv, x1, newton_param );
		f2 = eval_f_optimization ( load_f, v, dv, x2, newton_param );
	}

	if ( f_origin < f_optima )
	{
		x_optima = x_origin;
		if ( newton_param->debug )
		{
			printf( "[line_search] non-optimizion\n" );
		}
	}
	else
	{
		if ( newton_param->debug )
		{
			printf( "[line_search] reduce f_norm(1)=%.10le --> f_norm(%.15le)=%.10le\n", f_origin, x_optima, f_optima );
		}
	}


	return x_optima;
}

// <-----   L  ----->
// a____x₁____x₂____b
// 
// assume problem is convex hull, then
// x₁ = b - 𝜏*L
// x₂ = a + 𝜏*L
// if f(x₂) > f(x₁) then
// pruning x₂ ~ b and x₂ = a + 𝜏*(L - (b - (a + 𝜏*L))) = a + 𝜏*𝜏*L
// then solve algebra b - 𝜏*L = a + 𝜏*𝜏*L
// L*𝜏² + L*𝜏 - L = 𝜏² + 𝜏 - 1 = 0
// 𝜏 = (-1 ± √sqrt(5)) / 2 
// 𝜏 = (-1 + sqrt(5)) / 2 ≈ 0.618 (golden ratio)
// if f(x₁) < f(x₂) then prun x₂ = x₁, x₁ = b - 𝜏*L
// if f(x₂) < f(x₁) then prun x₁ = x₂, x₂ = a + 𝜏*L
// if f(x₁) = f(x₂) then prun a ~ x₁ and x₂ ~ b
// 
// almostly each iteration call only one f(x)
// linear converge rate fast then interval halving
//
static double golden_section ( void (load_f) (double *v, double*f), double *v, double *dv, double a, double b, newton_param_t *newton_param )
{
	int iter;
	double x_optima = NAN;
	double f_optima = NAN;
	double f1;
	double f2;
	double x1;
	double x2;
	double a_last;
	double b_last;
	double x_origin;
	double f_origin;
	double l;
	double delta;
	double golden_ratio = (-1 + sqrt(5)) / 2;
	double tol = newton_param->line_search_tol;

	x_origin = b;
	f_origin = eval_f_optimization ( load_f, v, dv, b, newton_param );

	iter = 1;
	l = b - a;
	delta = l * golden_ratio;
	x1 = b - delta;
	x2 = a + delta;
	f1 = eval_f_optimization ( load_f, v, dv, x1, newton_param );
	f2 = eval_f_optimization ( load_f, v, dv, x2, newton_param );

	if ( f1 <= f2 )
	{
		x_optima = x1;
		f_optima = f1;
	}
	else if ( f2 < f1 )
	{
		x_optima = x2;
		f_optima = f2;
	}

	while( l > tol ) 
	{
		if ( newton_param->debug )
		{
			printf( "[line_search] i=%d l=%.10le a=%.10le b=%.10le x1=%.10le x2=%.10le f1=%.10le f2=%.10le f_optima=%.10le\n", iter, l, a, b, x1, x2, f1, f2, f_optima );
		}

		if ( f1 < f2 )
		{
			if ( f1 < f_optima )
			{
				x_optima = x1;
				f_optima = f1;
			}
			b = x2;
			x2 = x1;
			f2 = f1;

			l = b - a;
			delta = l * golden_ratio;
			x1 = b - delta;
			f1 = eval_f_optimization ( load_f, v, dv, x1, newton_param );
		}
		else if ( f2 < f1 )
		{
			if ( f2 < f_optima )
			{
				x_optima = x2;
				f_optima = f2;
			}
			a = x1;
			x1 = x2;
			f1 = f2;

			l = b - a;
			delta = l * golden_ratio;
			x2 = a + delta;
			f2 = eval_f_optimization ( load_f, v, dv, x2, newton_param );
		}
		else
		{
			a = x1;
			b = x2;

			if ( f1 < f_optima )
			{
				x_optima = x1;
				f_optima = f1;
			}

			l = b - a;
			delta = l * golden_ratio;
			x1 = b - delta;
			x2 = a + delta;
			f1 = eval_f_optimization ( load_f, v, dv, x1, newton_param );
			f2 = eval_f_optimization ( load_f, v, dv, x2, newton_param );
		}

		++iter;

		l = b - a;
	}

	if ( f_origin < f_optima )
	{
		x_optima = x_origin;
		if ( newton_param->debug )
		{
			printf( "[line_search] non-optimizion\n" );
		}
	}
	else
	{
		if ( newton_param->debug )
		{
			printf( "[line_search] reduce f_norm(1)=%.10le --> f_norm(%.15le)=%.10le\n", f_origin, x_optima, f_optima );
		}
	}


	return x_optima;
}

static void dump_debug_data ( FILE *fout, int n, int iter, double *x, double *dx, double *f )
{
	if ( fout )
	{
		fprintf( fout, "%d ", iter );
		for ( int i = 0; i < n; ++i )
		{
			fprintf( fout, "%.15le ", x[i] );
		}
		for ( int i = 0; i < n; ++i )
		{
			fprintf( fout, "%.15le ", dx[i] );
		}
		for ( int i = 0; i < n; ++i )
		{
			fprintf( fout, "%.15le ", f[i] );
		}
		fprintf( fout, "\n" );
	}
}

static double eval_arc_length_constrain ( int n, double *x, double *xc, double p, double pc, double dt )
{
	double g = 0;
	for ( int i = 0; i < n; ++i )
	{
		g += (x[i] - xc[i]) * (x[i] - xc[i]);
	}
	g += (p - pc) * (p - pc);
	g -= dt * dt;
	return g;
}
