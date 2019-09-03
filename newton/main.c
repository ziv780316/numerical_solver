#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <dlfcn.h>

#include "newton.h"
#include "matrix_solver.h"
#include "opts.h"

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

		void (*p_load_f) (double *, double *) = (void (*)(double *, double *)) dlsym ( handle, "load_f" );
		if ( !p_load_f )
		{
			fprintf( stderr, "[Error] load symbol 'load_f' fail --> %s\n", dlerror() );
			abort();
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
		if ( NEWTON_CHORD_WITH_BYPASS_CHECK == g_opts.iterative_type )
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

		newton_iterative_type iterative_type = g_opts.iterative_type;
		newton_damped_type damped_type = g_opts.damped_type;
		newton_rescue_type rescue_type = g_opts.rescue_type;
		newton_derivative_type derivative_type = g_opts.diff_type;
		int maxiter = g_opts.maxiter;
		int miniter = g_opts.miniter;
		performance_stat nr_stat;
		memset( &nr_stat, 0, sizeof( performance_stat ) );
		double *x_result = (double *) malloc ( sizeof(double) * n );
		double *f_result = (double *) malloc ( sizeof(double) * n );
		double rtol = g_opts.rtol;
		double atol = g_opts.atol;
		double bypass_rtol = g_opts.bypass_rtol;
		double bypass_atol = g_opts.bypass_atol;
		double residual_tol = g_opts.residual_tol;
		double max_dx = g_opts.max_dx;
		double jmin = g_opts.jmin;
		bool random_initial = g_opts.random_initial;
		bool debug = g_opts.debug;
		bool converge;
		char *debug_file = g_opts.output_file;

		converge = newton_solve ( iterative_type, 
			       		  damped_type,
			       		  rescue_type,
			       		  derivative_type,
					  n,
			       		  x_init,
			       		  x_result,
			       		  f_result,
			       		  p_load_f,
			       		  p_load_jacobian,
					  p_bypass_check,
			       		  maxiter,
			       		  miniter,
			       		  rtol,
			       		  atol,
			       		  bypass_rtol,
			       		  bypass_atol,
			       		  residual_tol,
					  max_dx,
					  jmin,
			       		  random_initial,
					  &nr_stat,
			       		  debug,
					  debug_file );

		if ( !debug )
		{
			printf( "\n========== Newton Converge %s in %d Iteration ==========\n", (converge ? "Success" : "Fail"), nr_stat.n_iter );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.10e  f=%.10e\n", i, x_result[i], f_result[i] );
			}
		}
		printf( "* performance summary:\n" );
		printf( "n_iter       = %d\n", nr_stat.n_iter );
		printf( "n_mat_factor = %d\n", nr_stat.n_mat_factor );
		printf( "n_mat_solve  = %d\n", nr_stat.n_mat_solve );
		printf( "n_f_load     = %d\n", nr_stat.n_f_load );
		printf( "n_jac_load   = %d\n", nr_stat.n_jac_load );

		free( x_result );
		free( f_result );
		dlclose( handle );
	}

	return EXIT_SUCCESS;
}




