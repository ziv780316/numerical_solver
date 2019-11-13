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

		newton_param_t *newton_param = &(g_opts.newton_param);
		double *J;
		double *x_result = (double *) malloc ( sizeof(double) * n );
		double *f_result = (double *) malloc ( sizeof(double) * n );
		bool converge;
		char *debug_file = g_opts.output_file;

		converge = newton_solve ( newton_param,
			       		  J,
			       		  x_init,
			       		  px_ans,
			       		  x_result,
			       		  f_result,
			       		  p_load_f,
			       		  p_load_jacobian,
					  p_bypass_check,
					  debug_file );

		if ( !g_opts.newton_param.debug )
		{
			printf( "\n========== Newton Converge %s in %d Iteration ==========\n", (converge ? "Success" : "Fail"), newton_param->nr_stat.n_iter );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.10e  f=%.10e\n", i, x_result[i], f_result[i] );
			}
		}
		printf( "* performance summary:\n" );
		printf( "n_iter       = %d\n", newton_param->nr_stat.n_iter );
		printf( "n_mat_factor = %d\n", newton_param->nr_stat.n_mat_factor );
		printf( "n_mat_solve  = %d\n", newton_param->nr_stat.n_mat_solve );
		printf( "n_f_load     = %d\n", newton_param->nr_stat.n_f_load );
		printf( "n_jac_load   = %d\n", newton_param->nr_stat.n_jac_load );

		free( x_result );
		free( f_result );
		dlclose( handle );
	}

	return EXIT_SUCCESS;
}




