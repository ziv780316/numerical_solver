#include <stdio.h>
#include <stdlib.h>
#include "methods.h"
#include "opts.h"

opt_t g_opts;

int main( int argc, char **argv ) 
{
	// getopt parse command line options
	parse_cmd_options( argc, argv );

	bool debug = g_opts.debug;
	INTEGRATION_TYPE method = g_opts.method;
	int maxord = g_opts.maxord;
	int order = 1; // start-up order can only be 1
	double *ylist  = (double *) calloc ( maxord + 1, sizeof(double) );

	// initial solution
	double y0 = g_initial_solution;

	// time zone setting
	double h = g_opts.tstep;
	double xn_1 = 0;
	double yn_1 = y0;
	double yn = y0;
	double xstop = g_opts.tstop;
	
	// benchmark summary
	int total_points = 0;
	double lte;

	// imediately flush stream into termial
	setbuf( stdout, 0 );

	if ( debug )
	{
		printf( "x\ty\ty_exact\tlte\tlte(%%)\tdiff\n" );
	}
	while ( xn_1 <= xstop )
	{
		yn_1 = yn;

		if ( debug )
		{
			if ( 0 == total_points )
			{	
				printf( "%.10e %.10e %.10e %.10e %.10e %.10e\n", xn_1, yn_1, y0, 0.0, 0.0, diff(yn_1, xn_1) );
			}
			else
			{
				lte = exact(xn_1) - yn_1;
				printf( "%.10e %.10e %.10e %.10e %.10e %.10e\n", xn_1, yn_1, exact(xn_1), lte, 100.0 * (lte / exact(xn_1)), diff(yn_1, xn_1) );
			}
		}
		else
		{
			printf( "%.10e %.10e\n", xn_1, yn_1 );
		}

		// shift solution states
		for ( int i = maxord; i > 0; --i )
		{
			ylist[i] = ylist[i - 1];
		}

		switch ( method )
		{
			case AM:
			break;

			case AB:
			yn = adams_bashforth ( order, xn_1, yn_1, h, ylist );
			break;

			case BDF:
			break;

			case SIMPSON:
			// not implement yet...
			break;

			case RK:
			// not implement yet...
			break;
		}

		// go next time
		xn_1 += h;
		++total_points;
		if ( order < maxord )
		{
			++order;
		}
	}

	return 0;
}



