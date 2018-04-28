#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "methods.h"
#include "opts.h"

opt_t g_opts;
int g_total_iteration = 0;
int g_total_points = 0;
FILE *g_fout_local_solution = NULL;

int main( int argc, char **argv ) 
{
	// imediately flush stream into termial
	setbuf( stdout, 0 );

	// getopt parse command line options
	parse_cmd_options( argc, argv );

	bool debug = g_opts.debug;
	INTEGRATION_TYPE method = g_opts.method;
	int maxord = g_opts.maxord;
	int order = 1; // start-up order can only be 1
	double *ylist  = (double *) calloc ( maxord + 2, sizeof(double) ); // + 2 for polynomial extrapolation

	// initial solution
	double y0 = g_initial_solution;

	// time zone setting
	double h = g_opts.tstep;
	double tn_1 = 0;
	double yn_1 = y0;
	double yn = y0;
	double tstop = g_opts.tstop;
	
	// benchmark summary
	double lte;

	// linearize solution at each time point
	if ( g_opts.analysis_file )
	{
		g_fout_local_solution = fopen( g_opts.analysis_file, "w" );
		if ( !g_fout_local_solution )
		{
			fprintf( stderr, "[Error] open output file %s error --> %s\n", g_opts.analysis_file, strerror( errno ) );
			abort();
		}
		fprintf( g_fout_local_solution, "step(t,a) = (t >= a) ? 1 : 1/0\n" );
	}

	// start-up
	switch ( method )
	{
		case AM:
		case AB:
			ylist[0] = diff(y0, 0); 
			break;

		case BDF:
			ylist[0] = y0;
			break;

		case SIMPSON:
			// not implement yet...
			break;

		case RK:
			// not implement yet...
			break;
	}

	if ( debug )
	{
		printf( "t\ty\ty_exact\tlte\tlte(%%)\tdiff\n" );
	}
	while ( tn_1 <= tstop )
	{
		yn_1 = yn;

		// local analysis for eigenvalue estimation
		if ( g_opts.analysis_file )
		{
			double lamda = jacobian( yn_1, tn_1 );
			double q_over_lamda = diff( yn_1, tn_1 ) / lamda;

			fprintf( g_fout_local_solution, "t=%.10e lamda=%.10e\n", tn_1, lamda );
			fprintf( g_fout_local_solution, "f%d(t)=((%.10e)*exp(%.10e*(t - %.10e)) - %.10e + %.10e)*step(t,%.10e)\n", g_total_points, q_over_lamda, lamda, tn_1, q_over_lamda, yn_1, tn_1 );
		}

		if ( debug )
		{
			if ( 0 == g_total_points )
			{	
				printf( "%.10e %.10e %.10e %.10e %.10e %.10e\n", tn_1, yn_1, y0, 0.0, 0.0, diff(yn_1, tn_1) );
			}
			else
			{
				lte = exact(tn_1) - yn_1;
				printf( "%.10e %.10e %.10e %.10e %.10e %.10e\n", tn_1, yn_1, exact(tn_1), lte, 100.0 * (lte / exact(tn_1)), diff(yn_1, tn_1) );
			}
		}
		else
		{
			printf( "%.10e %.10e\n", tn_1, yn_1 );
		}

		// shift solution states, ylist[0] will be update in interation function
		for ( int i = maxord + 1; i > 0; --i )
		{
			ylist[i] = ylist[i - 1];
		}

		switch ( method )
		{
			case AM:
			yn = adams_moulton ( order, tn_1, yn_1, h, ylist );
			break;

			case AB:
			yn = adams_bashforth ( order, tn_1, yn_1, h, ylist );
			break;

			case BDF:
			yn = bdf ( order, tn_1, yn_1, h, ylist );
			break;

			case SIMPSON:
			// not implement yet...
			break;

			case RK:
			// not implement yet...
			break;
		}

		// show solution of y' of implicit method
		if ( g_opts.analysis_file )
		{
			if ( (BDF == method) && (1 == order) )
			{
				double tn = tn_1 + h;
				fprintf( g_fout_local_solution, "time:%.10e y'=%.10e, BDF%d y'=%.10e\n", tn, diff(yn, tn), order, (yn - yn_1)/h );
			}
		}

		// go next time
		tn_1 += h;
		++g_total_points;
		//if ( order < maxord && g_total_points > 1 ) // SPICE start up
		if ( order < maxord )
		{
			++order;
		}
	}

	if ( debug )
	{
		printf( "[Summary]\nTotal iterations = %d\nTotal timepoints = %d\n", g_total_iteration, g_total_points );
	}

	return 0;
}



