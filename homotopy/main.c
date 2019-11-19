#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <dlfcn.h>
#include <math.h>

#include "newton.h"
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

		double *px_ans = (double *) dlsym ( handle, "x_ans" );
		if ( !px_ans )
		{
			px_ans = NULL;
			dlerror(); 
		}

		double *plamda = (double *) dlsym ( handle, "lamda" );
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

		void (*p_load_df_dp) (double *, double *, double ) = (void (*)(double *, double *, double)) dlsym ( handle, "load_df_dp" );
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

		bool (*p_bypass_check) (double *, double *, double *) = NULL;
		if ( NEWTON_CHORD_WITH_BYPASS_CHECK == g_opts.newton_param.iterative_type )
		{
			p_bypass_check = (bool (*)(double *, double *, double *)) dlsym ( handle, "bypass_check" );
			if ( !p_bypass_check )
			{
				fprintf( stderr, "[Warning] load symbol 'bypass_check' fail --> %s\n", dlerror() );
				fprintf( stderr, "[Warning] use default __bypass_check routine\n" );
			}
			dlerror(); // clear error
		}

		// read initial_x0_file
		int n = *pn;
		double *x_init;
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
			x_init = px0;
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
		double *f_result = (double *) malloc ( sizeof(double) * n );
		double *f_delta = (double *) malloc ( sizeof(double) * n );
		double *df_dp = (double *) malloc ( sizeof(double) * n );
		double *dx_dp = (double *) malloc ( sizeof(double) * n );
		double *dx_dp_difference = (double *) malloc ( sizeof(double) * n );
		double *s0 = (double *) malloc ( sizeof(double) * n );
		double *s1 = (double *) malloc ( sizeof(double) * n );
		double *s_extrapolate = (double *) malloc ( sizeof(double) * n );
		double *s_extrapolate_difference = (double *) malloc ( sizeof(double) * n );
		double p0 = NAN;
		double p1 = NAN;
		double dp;
		bool converge;
		bool use_extrapolation = false;

		double lamda = 0;
		double lamda_old = 0;
		double lamda_start = 1e-3;	
		double lamda_end = 1;	

		while ( lamda <= lamda_end )
		{
			*plamda = lamda;
			++homotopy_param->hom_stat.n_step; 

			if ( homotopy_param->hom_stat.n_success >= 2)
			{
				if ( HOMOTOPY_EXTRAPOLATE_NONE == homotopy_param->extrapolate_type )
				{
					memcpy( x_init, s0, sizeof(double) * n );
				}
				else
				{
					use_extrapolation = true;
					dp = lamda - p0;
					for ( int i = 0; i < n; ++i )
					{
						dx_dp_difference[i] = (s0[i] - s1[i]) / (p0 - p1);
					}

					if ( HOMOTOPY_EXTRAPOLATE_DIFFERENCE == homotopy_param->extrapolate_type )
					{
						for ( int i = 0; i < n; ++i )
						{
							s_extrapolate[i] = s0[i] + dx_dp_difference[i] * dp;
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
							p_load_df_dp( s0, df_dp, p0 );
							++(homotopy_param->hom_stat.n_df_dp_load);
						}
						else
						{
							// forward difference to approximate ∂f/∂p
							double delta_ratio = 1e-6;
							double delta_lamda = lamda * delta_ratio;
							*plamda = p0 + delta_lamda;
							p_load_f( s0, f_delta );
							++(homotopy_param->hom_stat.n_f_load_sensitivity);
							*plamda = lamda;
							for ( int i = 0; i < n; ++i )
							{
								df_dp[i] = (f_delta[i] - f_result[i]) / delta_lamda;
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

						// xₖ₊₁ = xₖ + ∂x/∂p * Δp
						for ( int i = 0; i < n; ++i )
						{
							s_extrapolate[i] = s0[i] + dx_dp[i] * dp;
							s_extrapolate_difference[i] = s0[i] + dx_dp_difference[i] * dp;
						}
						if ( homotopy_param->debug )
						{
							printf( "* Sensitivity ∂x/∂λ\n" );
							for ( int i = 0; i < n; ++i )
							{
								printf( "%d: ∂x/∂λ=%.10e Δx/Δλ=%.10le\n", i, dx_dp[i], dx_dp_difference[i] );
							}
							printf( "\n" );
						}
					}

					memcpy( x_init, s_extrapolate, sizeof(double) * n );

					if ( homotopy_param->debug )
					{
						printf( "* Extrapolation p0=%.10le p1=%.10le\n", p0, p1 );
						for ( int i = 0; i < n; ++i )
						{
							printf( "x[%d]=%.10e  xpred[%d]=%.10le\n", i, s0[i], i, s_extrapolate[i] );
						}
						printf( "\n" );
					}
				}
			}

			memset( &(newton_param->nr_stat), 0, sizeof( performance_stat_t ) );

			converge = newton_solve ( newton_param,
					perm,
					J,
					x_init,
					px_ans,
					x_result,
					f_result,
					p_load_f,
					p_load_jacobian,
					p_bypass_check,
					NULL );

			update_homotopy_stat( homotopy_param, newton_param );

			printf( "\n* λ=%.10le NR Converge %s in %d Iteration\n", lamda, (converge ? "Success" : "Fail"), newton_param->nr_stat.n_iter );
			if ( homotopy_param->debug )
			{
				for ( int i = 0; i < n; ++i )
				{
					printf( "x[%d]=%.10e  f=%.10e\n", i, x_result[i], f_result[i] );
				}
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
								printf( "%d: sol=%.10e s0=%.10le xp=%.10le xpd=%.10le\n", i, x_result[i], s0[i], s_extrapolate[i], s_extrapolate_difference[i] );
								printf( "%d: dnone=%.10e dpred=%.10le ddpred=%.10le dratio=%.10g%% ddratio=%.10g%%\n", i, diff_none, diff_pred, diff_pred_difference, ratio_none - ratio_pred, ratio_none - ratio_pred_difference );
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
								printf( "%d: sol=%.10e s0=%.10le xp=%.10le\n", i, x_result[i], s0[i], s_extrapolate[i] );
								printf( "%d: diff_none=%.10e diff_pred=%.10le enhace_ratio=%.10g%%\n", i, diff_none, diff_pred, ratio_none - ratio_pred );
							}
						}
						printf( "\n" );
					}
				}

				++homotopy_param->hom_stat.n_success; 
				memcpy( s1, s0, sizeof(double) * n );
				memcpy( s0, x_result, sizeof(double) * n );
				p1 = p0;
				p0 = lamda;

				// ---------------------------------
				// use to output raw data for plot and analysis converge issue
				// ---------------------------------
				if ( fout_debug )
				{
					fprintf( fout_debug, "%.15le ", lamda );
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
				++homotopy_param->hom_stat.n_fail; 
			}

			if ( 0 == lamda )
			{
				if ( converge )
				{
					lamda = lamda_start;
				}
				else
				{
					printf( "homotopy fail in λ=0\n" );
					break;
				}
			}
			else
			{
				if ( converge )
				{
					if ( 1 == lamda )
					{
						// complete homotopy
						break;
					}
					else
					{
						lamda_old = lamda;
						lamda *= 1.1; // increase 10 %
						if ( lamda > 1 )
						{
							lamda = 1;
						}
					}
				}
				else
				{
					lamda = lamda * 0.55 + lamda_old * 0.45; // at lease lamda = lamda_old
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

		free( x_result );
		free( f_result );
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


