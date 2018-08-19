#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "newton.h"
#include "matrix_solver.h"
#include "opts.h"

int main ( int argc, char **argv )
{
	if ( 1 == argc )
	{
		newton_iterative_type iterative_type= NEWTON_NORMAL;
		newton_derivative_type derivative_type = NEWTON_DIFF_FORWARD;
		int maxiter = -1;
		double rtol = 1e-3;
		double atol = 1e-6;
		double residual_tol = 1e-9;
		bool random_initial = false;
		bool debug = true;
		newton_solve ( iterative_type, derivative_type, maxiter, rtol, atol, residual_tol, random_initial, debug );
	}
	else
	{
		parse_cmd_options ( argc, argv );
		newton_solve ( g_opts.iterative_type, g_opts.diff_type, g_opts.maxiter, g_opts.rtol, g_opts.atol, g_opts.residual_tol, g_opts.random_initial, g_opts.debug );
	}




	return EXIT_SUCCESS;
}




