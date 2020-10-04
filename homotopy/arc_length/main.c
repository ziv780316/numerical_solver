#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <dlfcn.h>
#include <math.h>

#include "matrix_solver.h"
#include "opts.h"

void update_homotopy_stat( homotopy_param_t *homotopy_param, newton_param_t *newton_param );

int main ( int argc, char **argv )
{
	if ( 1 == argc )
	{
		show_help();
	}
	else
	{
		// getopt parse command line arguments
		parse_cmd_options ( argc, argv );

		// open test function in *.so
		if ( !g_opts.problem_so )
		{
			fprintf( stderr, "[Error] please specify test function by '-p *.so'\n" );
			abort();
		}

		// dynamic link test function
		void *handle = dlopen( g_opts.problem_so, RTLD_NOW | RTLD_GLOBAL );
		if ( !handle )
		{
			fprintf( stderr, "[Error] dlopen '%s' fail --> %s\n", g_opts.problem_so, dlerror() );
			abort();
		}
		dlerror(); // clear error

		int *pn = (int *) dlsym ( handle, "nf" );
		if ( !pn )
		{
			fprintf( stderr, "[Error] load symbol 'nf' fail --> %s\n", dlerror() );
			abort();
		}

		double *px0 = (double *) dlsym ( handle, "x0" );
		if ( !px0 )
		{
			fprintf( stderr, "[Error] load symbol 'x0' fail --> %s\n", dlerror() );
			abort();
		}

		double *pp = (double *) dlsym ( handle, "p" );
		if ( !pp )
		{
			fprintf( stderr, "[Error] load symbol 'p' fail --> %s\n", dlerror() );
			abort();
		}

		double *px_ans = (double *) dlsym ( handle, "x_ans" );
		if ( !px_ans )
		{
			px_ans = NULL;
			dlerror(); 
		}

		double *plamda = (double *) dlsym ( handle, "p" );
		if ( !plamda )
		{
			plamda = NULL;
			dlerror(); 
		}

		void (*p_load_f) (double *, double *) = (void (*)(double *, double *)) dlsym ( handle, "load_f" );
		if ( !p_load_f )
		{
			fprintf( stderr, "[Error] load symbol 'load_f' fail --> %s\n", dlerror() );
			abort();
		}
		dlerror(); // clear error

		void (*p_load_df_dp) (double *, double *) = (void (*)(double *, double *)) dlsym ( handle, "load_df_dp" );
		if ( !p_load_df_dp )
		{
			fprintf( stderr, "[Warning] load symbol 'load_df_dp' fail --> %s, use forward difference approximate\n", dlerror() );
		}
		dlerror(); // clear error

		void (*p_load_jacobian) (double *, double *) = (void (*)(double *, double *)) dlsym ( handle, "load_jacobian" );
		if ( !p_load_jacobian )
		{
			fprintf( stderr, "[Error] load symbol 'load_jacobian' fail --> %s\n", dlerror() );
			abort();
		}
		dlerror(); // clear error

		// read initial_x0_file
		int n = *pn;
		double *x_init = (double *) malloc ( sizeof(double) * n );
		g_opts.newton_param.n = n;
		if ( g_opts.initial_x0_file )
		{
			FILE *fin = fopen( g_opts.initial_x0_file, "r" );	
			if ( !fin )
			{
				fprintf( stderr, "[Error] fopen '%s' fail --> %s\n", g_opts.initial_x0_file, strerror(errno) );
				abort();
			}

			int ret;
			x_init = (double *) malloc ( sizeof(double) * n );
			for ( int i = 0; i < n; ++i )
			{
				ret = fscanf( fin, "%lf", &(x_init[i]) );
				if ( EOF == ret && ((i + 1) <= n) )
				{
					fprintf( stderr, "[Error] initial x0 length %d < %d\n",  i, n );
					abort();
				}
			}
		}
		else
		{
			memcpy( x_init, px0, sizeof(double) * n );
		}

		// use to output raw data for plot and analysis converge issue
		char *debug_file = g_opts.output_file;
		FILE *fout_debug = NULL;
		if ( g_opts.homotopy_param.debug && debug_file )
		{
			char debug_file_name[BUFSIZ] = {0};
			sprintf( debug_file_name, "%s.trace", debug_file );
			fout_debug = fopen( debug_file_name, "w" );
			if ( !fout_debug )
			{
				fprintf( stderr, "[Error] open debug file %s fail --> %s\n", debug_file_name, strerror(errno) );
				abort();
			}
			fprintf( fout_debug, "# lamda " );
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "x_%d ", i );
			}
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "sp_%d ", i );
			}
			for ( int i = 0; i < n; ++i )
			{
				fprintf( fout_debug, "s0_%d ", i );
			}
			fprintf( fout_debug, "\n" );
		}

		newton_param_t *newton_param = &(g_opts.newton_param);
		homotopy_param_t *homotopy_param = &(g_opts.homotopy_param);
		int *perm = (int *) malloc ( sizeof(int) * n );
		double *J = (double *) malloc ( sizeof(double) * (n * n) );
		double *x_result = (double *) malloc ( sizeof(double) * n );
		double g;
		double *f_result = (double *) malloc ( sizeof(double) * n );
		double *f_delta = (double *) malloc ( sizeof(double) * n );
		double *df_dp = (double *) malloc ( sizeof(double) * n );
		double dp_dt;
		double dp_dt_difference;
		double *dx_dp = (double *) malloc ( sizeof(double) * n );
		double *dx_dt_difference = (double *) malloc ( sizeof(double) * n );
		double *dx_dt = (double *) malloc ( sizeof(double) * n );
		double *s0 = (double *) malloc ( sizeof(double) * n );
		double *s1 = (double *) malloc ( sizeof(double) * n );
		double *s_extrapolate = (double *) malloc ( sizeof(double) * n );
		double *s_extrapolate_difference = (double *) malloc ( sizeof(double) * n );
		double p_extrapolate;
		double p_init = 0;
		double p = p_init;
		double p0 = NAN;
		double p1 = NAN;
		double dt = homotopy_param->arc_length;
		double dt0 = homotopy_param->arc_length;
		double dp;
		bool converge;
		bool use_extrapolation = false;

		// solve initial p0 
		memset( &(newton_param->nr_stat), 0, sizeof( performance_stat_t ) );
		int maxiter_origin = newton_param->maxiter;
		newton_param->maxiter = 150;

		converge = arc_length_bbd_newton_solve ( 
				newton_param,
				perm, // permuation for matrix ordering
				J, // jacobian
				pp,
				x_init, // initial x
				NULL,
				p_init, // initial p
				x_result, // final x
				NULL, // final p
				f_result, // final f(x)
				NULL, // final g(x)
				NULL, // hyper circle central (xc,pc)
				NAN, // hyper circle central (xc,pc)
				0, // arc length
				p_load_f,
				p_load_df_dp,
				p_load_jacobian,
				debug_file );

		newton_param->maxiter = maxiter_origin;
		update_homotopy_stat( homotopy_param, newton_param );
		printf( "\n* λ=%.10le NR Converge %s in %d Iteration\n", p, (converge ? "Success" : "Fail"), newton_param->nr_stat.n_iter );
		if ( homotopy_param->debug )
		{
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.10le  f=%.10le\n", i, x_result[i], f_result[i] );
			}
			printf( "p=%.10le g=%.10le\n", p, g );
			printf( "\n" );
		}

		if ( converge )
		{
			memcpy( s0, x_result, sizeof(double) * n );
			memcpy( s1, x_result, sizeof(double) * n );
			p0 = p;
			p1 = p0;
		}
		else
		{
			fprintf( stderr, "homotopy fail due to p0=%.15le solve fail\n", p );
			exit(1);
		}

		while ( p <= 1 )
		{
			++homotopy_param->hom_stat.n_step; 

			if ( HOMOTOPY_EXTRAPOLATE_NONE == homotopy_param->extrapolate_type )
			{
				memcpy( x_init, s0, sizeof(double) * n );
				p_init = p0;
			}
			else
			{
				use_extrapolation = true;
				for ( int i = 0; i < n; ++i )
				{
					dx_dt_difference[i] = (s0[i] - s1[i]) / dt0;
				}
				dp_dt_difference = (p0 - p1) / dt0;

				if ( HOMOTOPY_EXTRAPOLATE_DIFFERENCE == homotopy_param->extrapolate_type )
				{
					for ( int i = 0; i < n; ++i )
					{
						s_extrapolate[i] = s0[i] + dx_dt_difference[i] * dt;
					}
				}
				else if ( HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL == homotopy_param->extrapolate_type )
				{
					// use sensitivity to extrapolate 
					// 0 = (∂f/∂p) * ∂x/∂p + ∂f/∂p
					// ∂x/∂p = (∂f/∂p)⁻¹ * -∂f/∂p
					// 
					if ( p_load_df_dp && (HOMOTOPY_DF_DP_EXACT == homotopy_param->df_dp_type) )
					{
						*pp = p0;
						p_load_df_dp( s0, df_dp );
						++(homotopy_param->hom_stat.n_df_dp_load);
					}
					else
					{
						// forward difference to approximate ∂f/∂p
						double delta_ratio = 1e-6;
						double delta_p;
						if ( p0 == 0 )
						{
							delta_p = 1e-6;
						}
						else
						{
							delta_p = p0 * delta_ratio;
						}
						*pp = p0 + delta_p;
						p_load_f( s0, f_delta );
						++(homotopy_param->hom_stat.n_f_load_sensitivity);
						*pp = p0;
						for ( int i = 0; i < n; ++i )
						{
							df_dp[i] = (f_delta[i] - f_result[i]) / delta_p;
						}
					}

					// prepare right hand side
					memcpy( dx_dp, df_dp, sizeof(double) * n );
					for ( int i = 0; i < n; ++i )
					{
						dx_dp[i] *= -1;
					}

					// solve ∂x/∂p
					bool matrix_solve_ok = dense_solve ( n, 1, J, dx_dp, perm, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, REAL_NUMBER );
					++(homotopy_param->hom_stat.n_mat_solve_sensitivity);
					if ( !matrix_solve_ok )
					{
						fprintf( stderr, "[Error] sensitivity LU solve fail\n" );
						abort();
					}

					// Δt is incremental arc length
					// xₖ₊₁ = xₖ + ∂x/∂t * Δt
					// pₖ₊₁ = pₖ + ∂p/∂t * Δt
					// --------------------------
					// sensitivity derive:
					// dx² + dp² = dt² 
					// ‖∂x/∂t‖² + ‖∂p/∂t‖² = 1
					// ‖∂x/∂p * ∂p/∂t‖² + ‖∂p/∂t‖² = 1
					// (‖∂x/∂p‖² * ‖∂p/∂t‖²) + ‖∂p/∂t‖² = 1
					// (1 + ‖∂x/∂p‖²) * ‖∂p/∂t‖² = 1
					// ‖∂p/∂t‖² = 1 / (1 + ‖∂x/∂p‖²)
					// ∂p/∂t = + 1 / sqrt(1 + ‖∂x/∂p‖²)
					// ∂x/∂t = ∂x/∂p * ∂p/∂t
					double s_norm_2 = 0;
					for ( int i = 0; i < n; ++i )
					{
						s_norm_2 += dx_dp[i] * dx_dp[i];
					}
					dp_dt = 1 / sqrt(1 + s_norm_2);
					for ( int i = 0; i < n; ++i )
					{
						dx_dt[i] = dx_dp[i] * dp_dt;
					}

					for ( int i = 0; i < n; ++i )
					{
						s_extrapolate[i] = s0[i] + dx_dt[i] * dt;
						s_extrapolate_difference[i] = s0[i] + dx_dt_difference[i] * dt;
					}
					p_extrapolate = p0 + dp_dt * dt;
					if ( homotopy_param->debug )
					{
						printf( "\n==================== Sensitivity ====================\n" );
						for ( int i = 0; i < n; ++i )
						{
							printf( "%d: ∂x/∂t=%.10le Δx/Δt=%.10le ∂x/∂p=%.10le ∂f/∂p=%.10le\n", i, dx_dt[i], dx_dt_difference[i], dx_dp[i], df_dp[i] );
						}
						printf( "%d: ∂p/∂t=%.10le Δp/Δt=%.10le\n", n, dp_dt, dp_dt_difference );
						printf( "=====================================================\n" );
					}
				}

				memcpy( x_init, s_extrapolate, sizeof(double) * n );
				p_init = p_extrapolate;
				

				if ( homotopy_param->debug )
				{
					printf( "\n==================== Extrapolation ====================\n" );
					for ( int i = 0; i < n; ++i )
					{
						printf( "x[%d]=%.10le  xpred[%d]=%.10le\n", i, s0[i], i, s_extrapolate[i] );
					}
					printf( "p=%.10le  ppred=%.10le\n", p0, p_extrapolate );
					printf( "=====================================================\n" );
				}
			}

			memset( &(newton_param->nr_stat), 0, sizeof( performance_stat_t ) );

			printf( "* solve NR λ=%.10le ...\n", p );
			converge = arc_length_bbd_newton_solve ( 
					newton_param,
					perm, // permuation for matrix ordering
					J, // jacobian
					pp,
					x_init, // initial x
					NULL,
					p_init, // initial p
					x_result, // final x
					&p, // final p
					f_result, // final f(x)
					&g, // final g(x)
					s0, // xc hyper circle central (xc,pc)
					p0, // pc hyper circle central (xc,pc)
					dt, // arc length
					p_load_f,
					p_load_df_dp,
					p_load_jacobian,
					debug_file );

			update_homotopy_stat( homotopy_param, newton_param );

			printf( "\n* λ=%.10le NR Converge %s in %d Iteration\n", p, (converge ? "Success" : "Fail"), newton_param->nr_stat.n_iter );
			if ( homotopy_param->debug )
			{
				for ( int i = 0; i < n; ++i )
				{
					printf( "x[%d]=%.10le  f=%.10le\n", i, x_result[i], f_result[i] );
				}
				printf( "p=%.10le g=%.10le\n", p, g );
				printf( "\n" );
			}

			if ( converge )
			{
				if ( homotopy_param->debug && use_extrapolation )
				{
					if ( HOMOTOPY_EXTRAPOLATE_NONE != homotopy_param->extrapolate_type )
					{
						double diff_none;
						double diff_pred;
						double diff_pred_difference;
						double ratio_none;
						double ratio_pred;
						double ratio_pred_difference;
						printf( "* Extrapolation Enhance Ratio\n" );
						if ( HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL == homotopy_param->extrapolate_type )
						{
							for ( int i = 0; i < n; ++i )
							{
								diff_none = x_result[i] - s0[i];
								diff_pred = x_result[i] - s_extrapolate[i];
								diff_pred_difference = x_result[i] - s_extrapolate_difference[i];
								ratio_none = fabs(diff_none / x_result[i]) * 100;
								ratio_pred = fabs(diff_pred / x_result[i]) * 100;
								ratio_pred_difference = fabs(diff_pred_difference / x_result[i]) * 100;
								printf( "%d: sol=%.10le s0=%.10le xp=%.10le xpd=%.10le\n", i, x_result[i], s0[i], s_extrapolate[i], s_extrapolate_difference[i] );
								printf( "%d: dnone=%.10le dpred=%.10le ddpred=%.10le dratio=%.10g%% ddratio=%.10g%%\n", i, diff_none, diff_pred, diff_pred_difference, ratio_none - ratio_pred, ratio_none - ratio_pred_difference );
							}
						}
						else
						{
							for ( int i = 0; i < n; ++i )
							{
								diff_none = x_result[i] - s0[i];
								diff_pred = x_result[i] - s_extrapolate[i];
								ratio_none = fabs(diff_none / x_result[i]) * 100;
								ratio_pred = fabs(diff_pred / x_result[i]) * 100;
								printf( "%d: sol=%.10le s0=%.10le xp=%.10le\n", i, x_result[i], s0[i], s_extrapolate[i] );
								printf( "%d: diff_none=%.10le diff_pred=%.10le enhace_ratio=%.10g%%\n", i, diff_none, diff_pred, ratio_none - ratio_pred );
							}
						}
						printf( "\n" );
					}
				}

				++homotopy_param->hom_stat.n_success; 
				memcpy( s1, s0, sizeof(double) * n );
				memcpy( s0, x_result, sizeof(double) * n );
				p1 = p0;
				p0 = p;
				dt *= 2;
				dt = fmin( dt, homotopy_param->arc_length );

				// ---------------------------------
				// use to output raw data for plot and analysis converge issue
				// ---------------------------------
				if ( fout_debug )
				{
					fprintf( fout_debug, "%.15le ", p );
					for ( int i = 0; i < n; ++i )
					{
						fprintf( fout_debug, "%.15le ", x_result[i] );
					}
					for ( int i = 0; i < n; ++i )
					{
						fprintf( fout_debug, "%.15le ", s_extrapolate[i] );
					}
					for ( int i = 0; i < n; ++i )
					{
						fprintf( fout_debug, "%.15le ", s0[i] );
					}
					fprintf( fout_debug, "\n" );
				}
			}
			else
			{
				p = p0;
				++homotopy_param->hom_stat.n_fail; 
				dt /= 2;
				if ( dt < 1e-10 )
				{
					printf( "homotopy fail since dt=%.15le < 1e-10\n", dt );
					break;
				}
			}
		}

		printf( "* Newton performance summary:\n" );
		printf( "n_iter           = %d\n", homotopy_param->hom_stat.n_iter );
		printf( "n_mat_factor     = %d\n", homotopy_param->hom_stat.n_mat_factor );
		printf( "n_mat_solve      = %d\n", homotopy_param->hom_stat.n_mat_solve );
		printf( "n_mat_solve_sens = %d\n", homotopy_param->hom_stat.n_mat_solve_sensitivity );
		printf( "toal_mat_solve   = %d\n", homotopy_param->hom_stat.n_mat_solve + homotopy_param->hom_stat.n_mat_solve_sensitivity );
		printf( "n_f_load         = %d\n", homotopy_param->hom_stat.n_f_load );
		printf( "n_f_load_sens    = %d\n", homotopy_param->hom_stat.n_f_load_sensitivity );
		printf( "total_f_load     = %d\n", homotopy_param->hom_stat.n_f_load + homotopy_param->hom_stat.n_f_load_sensitivity );
		printf( "n_df_dp_load     = %d\n", homotopy_param->hom_stat.n_df_dp_load );
		printf( "n_jac_load       = %d\n", homotopy_param->hom_stat.n_jac_load );
		printf( "\n" );
		printf( "* Homotopy performance summary:\n" );
		printf( "n_step        = %d\n", homotopy_param->hom_stat.n_step );
		printf( "n_success     = %d\n", homotopy_param->hom_stat.n_success );
		printf( "n_fail        = %d\n", homotopy_param->hom_stat.n_fail );
		printf( "n_limit_point = %d\n", homotopy_param->hom_stat.n_limit_point );

		dlclose( handle );
	}

	return EXIT_SUCCESS;
}

void update_homotopy_stat( homotopy_param_t *homotopy_param, newton_param_t *newton_param )
{
	homotopy_param->hom_stat.n_iter += newton_param->nr_stat.n_iter;
	homotopy_param->hom_stat.n_mat_factor += newton_param->nr_stat.n_mat_factor;
	homotopy_param->hom_stat.n_mat_solve += newton_param->nr_stat.n_mat_solve;
	homotopy_param->hom_stat.n_f_load += newton_param->nr_stat.n_f_load;
	homotopy_param->hom_stat.n_jac_load += newton_param->nr_stat.n_jac_load;
}


