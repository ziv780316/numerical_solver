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
		newton_modified_type modified_type = g_opts.modified_type;
		newton_derivative_type derivative_type = g_opts.diff_type;
		int maxiter = g_opts.maxiter;
		int miniter = g_opts.miniter;
		int total_iter = 0;
		double *x_result = (double *) malloc ( sizeof(double) * n );
		double *f_result = (double *) malloc ( sizeof(double) * n );
		double rtol = g_opts.rtol;
		double atol = g_opts.atol;
		double residual_tol = g_opts.residual_tol;
		bool random_initial = g_opts.random_initial;
		bool debug = g_opts.debug;
		bool converge;

		converge = newton_solve ( iterative_type, 
			       		  modified_type,
			       		  derivative_type,
					  n,
			       		  x_init,
			       		  x_result,
			       		  f_result,
			       		  p_load_f,
			       		  p_load_jacobian,
			       		  maxiter,
			       		  miniter,
					  &total_iter,
			       		  rtol,
			       		  atol,
			       		  residual_tol,
			       		  random_initial,
			       		  debug );

		if ( !debug )
		{
			printf( "\n========== Newton Converge %s in %d Iteration ==========\n", (converge ? "Success" : "Fail"), total_iter );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.10e  f=%.10e\n", i, x_result[i], f_result[i] );
			}
		}

		free( x_result );
		free( f_result );
		dlclose( handle );
	}

	return EXIT_SUCCESS;
}




