// test homotopy arc-length method to trace equilibrium path H(x,λ)=0
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
			exit(1);
		}

		// dynamic link test function
		void *handle = dlopen( g_opts.problem_so, RTLD_NOW | RTLD_GLOBAL );
		if ( !handle )
		{
			fprintf( stderr, "[Error] dlopen '%s' fail --> %s\n", g_opts.problem_so, dlerror() );
			exit(1);
		}
		dlerror(); // clear error

		int *pn = (int *) dlsym ( handle, "nf" );
		if ( !pn )
		{
			fprintf( stderr, "[Error] load symbol 'nf' fail --> %s\n", dlerror() );
			exit(1);
		}

		double *px0 = (double *) dlsym ( handle, "x0" );
		if ( !px0 )
		{
			fprintf( stderr, "[Error] load symbol 'x0' fail --> %s\n", dlerror() );
			exit(1);
		}

		double *pp = (double *) dlsym ( handle, "p" );
		if ( !pp )
		{
			fprintf( stderr, "[Error] load symbol 'p' fail --> %s\n", dlerror() );
			exit(1);
		}

		double *pp0 = (double *) dlsym ( handle, "p0" );
		if ( !pp0 )
		{
			pp0 = NULL;
			dlerror(); 
		}

		double *px_ans = (double *) dlsym ( handle, "x_ans" );
		if ( !px_ans )
		{
			px_ans = NULL;
			dlerror(); 
		}

		void (*p_load_f) (double *, double *) = (void (*)(double *, double *)) dlsym ( handle, "load_f" );
		if ( !p_load_f )
		{
			fprintf( stderr, "[Error] load symbol 'load_f' fail --> %s\n", dlerror() );
			exit(1);
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
			exit(1);
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
				exit(1);
			}

			int ret;
			x_init = (double *) malloc ( sizeof(double) * n );
			for ( int i = 0; i < n; ++i )
			{
				ret = fscanf( fin, "%lf", &(x_init[i]) );
				if ( EOF == ret && ((i + 1) <= n) )
				{
					fprintf( stderr, "[Error] initial x0 length %d < %d\n",  i, n );
					exit(1);
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
		FILE *fout_debug_tangent = NULL;
		if ( g_opts.homotopy_param.debug && debug_file )
		{
			char debug_file_name[BUFSIZ] = {0};
			sprintf( debug_file_name, "%s.trace", debug_file );
			fout_debug = fopen( debug_file_name, "w" );
			if ( !fout_debug )
			{
				fprintf( stderr, "[Error] open debug file %s fail --> %s\n", debug_file_name, strerror(errno) );
				exit(1);
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
			fprintf( fout_debug, "is_backtrace" );
			fprintf( fout_debug, "\n" );

			sprintf( debug_file_name, "%s.tangent", debug_file );
			fout_debug_tangent = fopen( debug_file_name, "w" );
			if ( !fout_debug_tangent )
			{
				fprintf( stderr, "[Error] open debug file %s fail --> %s\n", debug_file_name, strerror(errno) );
				exit(1);
			}
		}

		newton_param_t *newton_param = &(g_opts.newton_param);
		homotopy_param_t *homotopy_param = &(g_opts.homotopy_param);
		int *perm = (int *) malloc ( sizeof(int) * n );
		int J_size = n*n;
		double *J = (double *) malloc ( sizeof(double) * J_size );
		double *A = (double *) malloc ( sizeof(double) * ((n+1) * (n+1)) );
		double *x_result = (double *) malloc ( sizeof(double) * n );
		double g;
		double det_a = NAN;
		double *f_result = (double *) malloc ( sizeof(double) * n );
		double *f_delta = (double *) malloc ( sizeof(double) * n );
		double *df_dp = (double *) malloc ( sizeof(double) * n );
		double dp_dt = NAN;
		double dp_dt_last = NAN;
		double dp_dt_difference = NAN;
		double *dx_dp = (double *) malloc ( sizeof(double) * n );
		double *dx_dp_0 = (double *) malloc ( sizeof(double) * n );
		double *dx_dt_difference = (double *) malloc ( sizeof(double) * n );
		double *dx_dt = (double *) malloc ( sizeof(double) * n );
		double *s0 = (double *) malloc ( sizeof(double) * n );
		double *s1 = (double *) malloc ( sizeof(double) * n );
		double *s_extrapolate = (double *) malloc ( sizeof(double) * n );
		double *s_extrapolate_difference = (double *) malloc ( sizeof(double) * n );
		double p_extrapolate = NAN;
		double p_init = (pp0 ? *pp0 : 0);
		double p = p_init;
		double p0 = NAN;
		double p1 = NAN;
		double dt = homotopy_param->arc_length;
		double dt0 = homotopy_param->arc_length;
		double dp;
		bool final_sim = false;
		bool converge;
		bool use_extrapolation = false;
		bool matrix_solve_ok = false;
		bool need_reverse_sign = false;
		bool is_backtrace = false;
		double forward_step_cross_product = NAN;
		bool *var_backtrace_status = (bool *) malloc ( sizeof(bool) * n );
		double *var_forward_det = (double *) malloc ( sizeof(double) * J_size );
		double *var_det = (double *) malloc ( sizeof(double) * J_size );
		double Asub[2][2];
		int *var_det_violate_cnt = (int *) malloc ( sizeof(int) * n );
		int total_var_det_violate_cnt;
		double violate_thr;
		double violate_ratio;

		// solve initial p0 
		memset( &(newton_param->nr_stat), 0, sizeof( performance_stat_t ) );
		int maxiter_origin = newton_param->maxiter;
		newton_param->maxiter = 150;

		converge = arc_length_bbd_newton_solve ( 
				homotopy_param,
				newton_param,
				perm, // permuation for matrix ordering
				J, // jacobian
				pp,
				x_init, // initial x
				NULL,
				p_init, // initial p
				NULL,
				NAN,
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

			if ( fout_debug )
			{
				fprintf( fout_debug, "0 " );
				for ( int i = 0; i < n; ++i )
				{
					fprintf( fout_debug, "%.15le ", x_result[i] );
				}
				for ( int i = 0; i < n; ++i )
				{
					fprintf( fout_debug, "0 " );
				}
				for ( int i = 0; i < n; ++i )
				{
					fprintf( fout_debug, "0 " );
				}
				fprintf( fout_debug, "%d", is_backtrace );
				fprintf( fout_debug, "\n" );
			}
		}
		else
		{
			fprintf( stderr, "homotopy fail due to p0=%.15le solve fail\n", p );
			exit(1);
		}

		while ( (p <= 1) && (homotopy_param->hom_stat.n_step < homotopy_param->maxsteps) )
		{
			++homotopy_param->hom_stat.n_step; 

			if ( p < 0 )
			{
				break;
			}

			if ( HOMOTOPY_EXTRAPOLATE_NONE == homotopy_param->extrapolate_type )
			{
				memcpy( x_init, s0, sizeof(double) * n );
				p_init = p0;
			}
			else
			{
				if ( converge )
				{
					use_extrapolation = true;
					for ( int i = 0; i < n; ++i )
					{
						dx_dt_difference[i] = (s0[i] - s1[i]) / dt0;
					}
					dp_dt_difference = (p0 - p1) / dt0;

					// do extrapolation
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
						memcpy( dx_dp_0, dx_dp, sizeof(double) * n );
						memcpy( dx_dp, df_dp, sizeof(double) * n );
						for ( int i = 0; i < n; ++i )
						{
							dx_dp[i] *= -1;
						}

						// solve ∂x/∂p
						matrix_solve_ok = dense_solve ( n, 1, J, dx_dp, perm, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, REAL_NUMBER );
						++(homotopy_param->hom_stat.n_mat_solve_sensitivity);
						if ( !matrix_solve_ok )
						{
							fprintf( stderr, "[Error] sensitivity LU solve fail\n" );
							exit(1);
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
						dp_dt = 1 / sqrt(1 + s_norm_2); // XXX mustbe positive, it's necessary to check sign later 

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
								printf( "%d: ∂x/∂t=%.10le Δx/Δt=%.10le ∂x/∂p=%.10le ∂f/∂p=%.10le tangent_line(p)=%.10le*(p-%.10le)+%.10le\n", i, dx_dt[i], dx_dt_difference[i], dx_dp[i], df_dp[i] , dx_dt[i]/dp_dt, p0, s0[i]);
							}
							printf( "%d: ∂p/∂t=%.10le Δp/Δt=%.10le\n", n, dp_dt, dp_dt_difference );
							printf( "=====================================================\n" );
						}
					}

					// backtrace handling
					if ( 0 == homotopy_param->hom_stat.n_success )
					{
						// do nothing
					}
					else if ( HOMOTOPY_ARC_LENGTH_BACKTRACE_HANDLE_CROSS_PRODUCT == homotopy_param->arc_length_backtrace_type )
					{
						printf( "\n* evaluate λ=%.10le cross product sign ...\n", p0 );

						// A₁₁ = ∂f/∂x = J(x)
						memset( J, 0, sizeof(double) * J_size );
						p_load_jacobian( s0, J );
						for ( int col = 0; col < n; ++col )
						{
							for ( int row = 0; row < n; ++row )
							{
								*(A + col*(n+1) + row) = *(J + col*n + row);
							}
						}

						// A₁₂ = ∂f/∂p
						p_load_df_dp( s0, df_dp );
						for ( int row = 0; row < n; ++row )
						{
							*(A + (n)*(n+1) + row) = df_dp[row];
						}

						// A₂₁ = ∂x/∂t
						if ( HOMOTOPY_EXTRAPOLATE_DIFFERENCE == homotopy_param->extrapolate_type )
						{
							for ( int col = 0; col < n; ++col )
							{
								*(A + col*(n+1) + n) = dx_dt_difference[col];
							}
						}
						else if ( HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL == homotopy_param->extrapolate_type )
						{
							for ( int col = 0; col < n; ++col )
							{
								*(A + col*(n+1) + n) = dx_dt[col];
							}
						}

						// A₂₂ = ∂p/∂t
						if ( HOMOTOPY_EXTRAPOLATE_DIFFERENCE == homotopy_param->extrapolate_type )
						{
							*(A + n*(n+1) + n) = dp_dt_difference;
						}
						else if ( HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL == homotopy_param->extrapolate_type )
						{
							*(A + n*(n+1) + n) = dp_dt;
						}


						if ( homotopy_param->debug )
						{
							printf( "A = \n" );
							dense_print_matrix( n + 1, n + 1, A, REAL_NUMBER );
							printf( "|A|=%.10le\n", det_a );
							det_a = dense_eval_det( n + 1, A, FACTOR_LU_RIGHT_LOOKING, REAL_NUMBER );
						}

						if ( 1 == homotopy_param->hom_stat.n_success )
						{
							forward_step_cross_product = det_a;
						}
						else
						{
							if ( det_a * forward_step_cross_product < 0 )
							{
								if ( HOMOTOPY_EXTRAPOLATE_DIFFERENCE == homotopy_param->extrapolate_type )
								{
								}
								else if ( HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL == homotopy_param->extrapolate_type )
								{
									dp_dt *= -1; // change tagent line sign
									for ( int i = 0; i < n; ++i )
									{
										dx_dt[i] = dx_dp[i] * dp_dt;
									}

									if ( homotopy_param->debug )
									{
										printf( "\n==================== Modified Sensitivity ====================\n" );
										for ( int i = 0; i < n; ++i )
										{
											printf( "%d: ∂x/∂t=%.10le Δx/Δt=%.10le ∂x/∂p=%.10le ∂f/∂p=%.10le tangent_line(p)=%.10le*(p-%.10le)+%.10le\n", i, dx_dt[i], dx_dt_difference[i], dx_dp[i], df_dp[i] , dx_dt[i]/dp_dt, p0, s0[i]);
										}
										printf( "%d: ∂p/∂t=%.10le Δp/Δt=%.10le\n", n, dp_dt, dp_dt_difference );
										printf( "=====================================================\n" );
										printf( "\n==================== Modified Extrapolation ====================\n" );
										for ( int i = 0; i < n; ++i )
										{
											printf( "xpred[%d]=%.10le -> %.10le\n", i, s_extrapolate[i], s0[i] + dx_dt[i] * dt );
										}
										printf( "ppred=%.10le -> %.10le\n", p_extrapolate, p0 + dp_dt * dt );
										printf( "=====================================================\n" );
									}
									for ( int i = 0; i < n; ++i )
									{
										s_extrapolate[i] = s0[i] + dx_dt[i] * dt;
										s_extrapolate_difference[i] = s0[i] + dx_dt_difference[i] * dt;
									}
									p_extrapolate = p0 + dp_dt * dt;
								}
							}
						}
					}
					else if ( (HOMOTOPY_ARC_LENGTH_BACKTRACE_HANDLE_SUB_CROSS_PRODUCT == homotopy_param->arc_length_backtrace_type) || 
						  (HOMOTOPY_ARC_LENGTH_BACKTRACE_HANDLE_DIAG_CROSS_PRODUCT == homotopy_param->arc_length_backtrace_type) )
					{
						memset( var_backtrace_status, 0, sizeof(bool) * n );
						memset( var_det_violate_cnt, 0, sizeof(int) * n );
						total_var_det_violate_cnt = 0;

						printf( "\n* evaluate λ=%.10le sub cross product sign ...\n", p0 );

						// Asub₂₁ = ∂x/∂t
						// Asub₂₂ = ∂p/∂t
						// find |Asub| for each variable and F
						p_load_jacobian( s0, J );
						p_load_df_dp( s0, df_dp );
						for ( int col = 0; col < n; ++col )
						{
							for ( int row = 0; row < n; ++row )
							{
						  		if ( HOMOTOPY_ARC_LENGTH_BACKTRACE_HANDLE_DIAG_CROSS_PRODUCT == homotopy_param->arc_length_backtrace_type )
								{
									if ( col != row )
									{
										continue;
									}
								}

								// Asub₁₁ = ∂f/∂x 
								Asub[0][0] = *(J + col*n + row);

								// Asub₁₂ = ∂f/∂p
								Asub[1][0] = df_dp[row];

								// Asub₂₁ = ∂x/∂t
								if ( HOMOTOPY_EXTRAPOLATE_DIFFERENCE == homotopy_param->extrapolate_type )
								{
									Asub[0][1] = dx_dt_difference[col];
								}
								else if ( HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL == homotopy_param->extrapolate_type )
								{
									Asub[0][1] = dx_dt[col];
								}
								else
								{
									Asub[0][1] = NAN;
								}

								// Asub₂₂ = ∂p/∂t
								if ( HOMOTOPY_EXTRAPOLATE_DIFFERENCE == homotopy_param->extrapolate_type )
								{
									Asub[1][1] = dp_dt_difference;
								}
								else if ( HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL == homotopy_param->extrapolate_type )
								{
									Asub[1][1] = dp_dt;
								}
								else
								{
									Asub[1][1] = NAN;
								}

								*(var_det + col*n + row) = (Asub[0][0]*Asub[1][1]) - (Asub[1][0]*Asub[0][1]);

								if ( homotopy_param->debug )
								{
									printf( "A%d%d = \n", row, col );
									printf( "%.10le %.10le\n", Asub[0][0], Asub[1][0] );
									printf( "%.10le %.10le\n", Asub[0][1], Asub[1][1] );
									printf( "|A%d%d|=%.10le\n", row, col, *(var_det + col*n + row) );
								}
							}
						}

						if ( 1 == homotopy_param->hom_stat.n_success )
						{
							memcpy( var_forward_det, var_det, sizeof(double) * J_size );
						}
						else
						{
							for ( int col = 0; col < n; ++col )
							{
								for ( int row = 0; row < n; ++row )
								{
									if ( HOMOTOPY_ARC_LENGTH_BACKTRACE_HANDLE_DIAG_CROSS_PRODUCT == homotopy_param->arc_length_backtrace_type )
									{
										if ( col != row )
										{
											continue;
										}
									}
									if ( *(var_det + col*n + row) * *(var_forward_det + col*n + row) < 0 )
									{
										var_det_violate_cnt[col] += 1;
										total_var_det_violate_cnt += 1;
									}
								}
							}

							violate_thr = 0.5;
							if ( HOMOTOPY_ARC_LENGTH_BACKTRACE_HANDLE_DIAG_CROSS_PRODUCT == homotopy_param->arc_length_backtrace_type )
							{
								violate_ratio = total_var_det_violate_cnt/(double)n;
							}
							else
							{
								violate_ratio = total_var_det_violate_cnt/(double)J_size;
							}
							if ( homotopy_param->debug )
							{
								if ( HOMOTOPY_ARC_LENGTH_BACKTRACE_HANDLE_DIAG_CROSS_PRODUCT == homotopy_param->arc_length_backtrace_type )
								{
									printf( "total %d/%d (%.3lf%%) |Asub| violation, thr=%.3lf%%\n", total_var_det_violate_cnt, n, violate_ratio, violate_thr );
								}
								else
								{
									printf( "total %d/%d (%.3lf%%) |Asub| violation, thr=%.3lf%%\n", total_var_det_violate_cnt, n*n, violate_ratio, violate_thr );
								}
							}

							if ( violate_ratio >= violate_thr )
							{
								dp_dt *= -1; // change tagent line sign
								for ( int i = 0; i < n; ++i )
								{
									dx_dt[i] = dx_dp[i] * dp_dt;
								}

								if ( homotopy_param->debug )
								{
									printf( "\n==================== Modified Sensitivity ====================\n" );
									for ( int i = 0; i < n; ++i )
									{
										printf( "%d: ∂x/∂t=%.10le Δx/Δt=%.10le ∂x/∂p=%.10le ∂f/∂p=%.10le tangent_line(p)=%.10le*(p-%.10le)+%.10le\n", i, dx_dt[i], dx_dt_difference[i], dx_dp[i], df_dp[i] , dx_dt[i]/dp_dt, p0, s0[i]);
									}
									printf( "%d: ∂p/∂t=%.10le Δp/Δt=%.10le\n", n, dp_dt, dp_dt_difference );
									printf( "=====================================================\n" );
									printf( "\n==================== Modified Extrapolation ====================\n" );
									for ( int i = 0; i < n; ++i )
									{
										printf( "xpred[%d]=%.10le -> %.10le\n", i, s_extrapolate[i], s0[i] + dx_dt[i] * dt );
									}
									printf( "ppred=%.10le -> %.10le\n", p_extrapolate, p0 + dp_dt * dt );
									printf( "=====================================================\n" );
								}
								for ( int i = 0; i < n; ++i )
								{
									s_extrapolate[i] = s0[i] + dx_dt[i] * dt;
									s_extrapolate_difference[i] = s0[i] + dx_dt_difference[i] * dt;
								}
								p_extrapolate = p0 + dp_dt * dt;
							}
						}
					}
					else if ( HOMOTOPY_ARC_LENGTH_BACKTRACE_HANDLE_DIFFERENCE == homotopy_param->arc_length_backtrace_type )
					{
						// turnning point reverse sign
						if ( dp_dt_difference < 0 )
						{
							if ( dp_dt > 0 )
							{
								printf( "* Turnning point is occured arround λ=%.10le\n", p1 );
								dp_dt *= -1;

								for ( int i = 0; i < n; ++i )
								{
									dx_dt[i] = dx_dp[i] * dp_dt;
								}
								if ( homotopy_param->debug )
								{
									printf( "\n==================== Modified Sensitivity ====================\n" );
									for ( int i = 0; i < n; ++i )
									{
										printf( "%d: ∂x/∂t=%.10le Δx/Δt=%.10le ∂x/∂p=%.10le ∂f/∂p=%.10le tangent_line(p)=%.10le*(p-%.10le)+%.10le\n", i, dx_dt[i], dx_dt_difference[i], dx_dp[i], df_dp[i] , dx_dt[i]/dp_dt, p0, s0[i]);
									}
									printf( "%d: ∂p/∂t=%.10le Δp/Δt=%.10le\n", n, dp_dt, dp_dt_difference );
									printf( "=====================================================\n" );
									printf( "\n==================== Modified Extrapolation ====================\n" );
									for ( int i = 0; i < n; ++i )
									{
										printf( "xpred[%d]=%.10le -> %.10le\n", i, s_extrapolate[i], s0[i] + dx_dt[i] * dt );
									}
									printf( "ppred=%.10le -> %.10le\n", p_extrapolate, p0 + dp_dt * dt );
									printf( "=====================================================\n" );
								}
								for ( int i = 0; i < n; ++i )
								{
									s_extrapolate[i] = s0[i] + dx_dt[i] * dt;
									s_extrapolate_difference[i] = s0[i] + dx_dt_difference[i] * dt;
								}
								p_extrapolate = p0 + dp_dt * dt;
							}
						}
					}
					else
					{
						// none
					}
				}
				else
				{
					if ( homotopy_param->debug )
					{
						printf( "\n==================== Sensitivity ====================\n" );
						for ( int i = 0; i < n; ++i )
						{
							printf( "%d: ∂x/∂t=%.10le Δx/Δt=%.10le ∂x/∂p=%.10le ∂f/∂p=%.10le tangent_line(p)=%.10le*(p-%.10le)+%.10le\n", i, dx_dt[i], dx_dt_difference[i], dx_dp[i], df_dp[i] , dx_dt[i]/dp_dt, p0, s0[i]);
						}
						printf( "%d: ∂p/∂t=%.10le Δp/Δt=%.10le\n", n, dp_dt, dp_dt_difference );
						printf( "=====================================================\n" );
					}

					// use last successful sensitivity results
					for ( int i = 0; i < n; ++i )
					{
						s_extrapolate[i] = s0[i] + dx_dt[i] * dt;
						s_extrapolate_difference[i] = s0[i] + dx_dt_difference[i] * dt;
					}
					p_extrapolate = p0 + dp_dt * dt;
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

			if ( HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL == homotopy_param->extrapolate_type )
			{
				if ( dp_dt_last * dp_dt < 0 )
				{
					printf( "* Turnning point is occured arround λ=%.10le\n", p1 );
				}
				dp_dt_last = dp_dt;
			}

			memset( &(newton_param->nr_stat), 0, sizeof( performance_stat_t ) );

			if ( p_init > 1 )
			{
				p_init = 1;
				p = p_init;
				final_sim = true;
				printf( "\n* Solve final λ=1 NR (step=%d) )...\n", homotopy_param->hom_stat.n_step );
				converge = arc_length_bbd_newton_solve ( 
						homotopy_param,
						newton_param,
						perm, // permuation for matrix ordering
						J, // jacobian
						pp,
						x_init, // initial x
						NULL,
						1, // final p
						NULL,
						NAN,
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
			}
			else
			{
				final_sim = false;
				p = p_init;
				is_backtrace = (p < p0);
				printf( "* Solve NR λ=%.10le (step=%d backtrace=%s) ...\n", p_init, homotopy_param->hom_stat.n_step, (is_backtrace? "yes" : "no") );
				converge = arc_length_bbd_newton_solve ( 
						homotopy_param,
						newton_param,
						perm, // permuation for matrix ordering
						J, // jacobian
						pp,
						x_init, // initial x
						NULL,
						p_init, // initial p
						dx_dt,
						dp_dt,
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
			}

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
				if ( homotopy_param->debug && fout_debug_tangent )
				{
					printf( "\n==================== Tangent Line ====================\n" );
					for ( int i = 0; i < n; ++i )
					{
						fprintf( fout_debug_tangent, "tangent_line_%d_%d(p)=%.10le*(p-%.10le)+%.10le\n", i, homotopy_param->hom_stat.n_success, dx_dt[i]/dp_dt, p0, s0[i]);
						fprintf( fout_debug_tangent, "delta_line_%d_%d(p)=%.10le*(p-%.10le)+%.10le\n", i, homotopy_param->hom_stat.n_success, (x_result[i]-s0[i])/(p - p0), p0, s0[i]);
					}
					printf( "=====================================================\n" );
				}

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
				dt = dt * 2;
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
					fprintf( fout_debug, "%d", is_backtrace );
					fprintf( fout_debug, "\n" );
				}

				if ( final_sim )
				{
					break;
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

		if ( final_sim && converge )
		{
			printf( "\n* Homotopy Success on λ=1 with total %d iterations\n", homotopy_param->hom_stat.n_iter );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.10le  f=%.10le\n", i, x_result[i], f_result[i] );
			}
		}
		else if ( p < 0 )
		{
			printf( "\n* Final λ=%.10le < 0, Homotopy Fail (there is no path to λ=1)\n", p );
		}
		else
		{
			printf( "\n* Homotopy Fail stuck on λ=%.10le with toal %d Iteration\n", p, homotopy_param->hom_stat.n_iter );
		}

		printf( "\n* Newton performance summary:\n" );
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


