#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "opts.h"
#include "methods.h"

double bdf ( int order, double tn_1, double yn_1, double hn_1, double *ylist )
{
	double yn;
	double yn_2 = ylist[2];
	double yn_3 = ylist[3];
	double tn = tn_1 + hn_1;
	double yn_predictor;
	double yn_k;
	double yn_k_delta;
	double diff_0;
	double diff_0_delta;

	// Newton-Raphoson method
	int iterno = 0;
	bool converge = false;
	double reltol = 0;
	double abstol = 1e-6;
	double nr_diff;
	double yn_next;
	double g;
	double g_delta;
	double g_diff;
	double delta = 1e-6;
	double *ylist_predictor = (double *) calloc ( order + 1, sizeof(double) );
	memcpy( ylist_predictor, ylist, sizeof(double) * (order + 1) );

	if ( g_opts.use_predictor )
	{
		// not implement yet
		yn = yn_1;
	}
	else
	{
		yn = yn_1; // NR initial value without predictor
	}

	do 
	{
		diff_0       = diff(yn, tn);
		diff_0_delta = diff(yn + delta, tn); // forward difference

		switch ( order )
		{
			case 1:
				// Backward-Euler
				yn_k = yn_1 + hn_1 * diff_0;
				yn_k_delta = yn_1 + hn_1 * diff_0_delta;
				break;
			case 2:
				// GEAR-2
				yn_k = (4.0 / 3.0) * yn_1 - (1.0 / 3.0) * yn_2 + (2.0 / 3.0) * hn_1 * diff_0;
				yn_k_delta = (4.0 / 3.0) * yn_1 - (1.0 / 3.0) * yn_2 + (2.0 / 3.0) * hn_1 * diff_0_delta;
				break;
			case 3:
				yn_k = (18.0 / 11.0) * yn_1 - (9.0 / 11.0) * yn_2 + (2.0 / 11.0) * yn_3 + (6.0 / 11.0) * hn_1 * diff_0;
				yn_k_delta = (18.0 / 11.0) * yn_1 - (9.0 / 11.0) * yn_2 + (2.0 / 11.0) * yn_3 + (6.0 / 11.0) * hn_1 * diff_0_delta;
				break;
			default:
				fprintf( stderr, "[Error] cannot support BDF order=%d\n", order );
				abort();
				break;
		}

		g = yn - yn_k;
		g_delta = (yn + delta) - yn_k_delta;
		g_diff = (g_delta - g) / delta; // forward difference
		yn_next = yn - (g / g_diff);
		nr_diff = fabs(yn_next - yn);
		yn = yn_next;
		++iterno;
		//printf( "%d: yn_next=%.10e yn=%.10e yn_1=%.10e yn_k=%.10e yn_predictor=%.10e diff_0=%.10e g=%.10e, g_delta=%.10e, g_diff=%.10e\n", iterno, yn_next, yn, yn_1, yn_k, yn_predictor, diff_0, g, g_delta, g_diff );

		if ( nr_diff > (reltol * fmax(fabs(yn), fabs(yn_next))) + abstol )
		{
			converge = false;
		}
		else
		{
			converge = true;
		}
	} while ( !converge );

	free( ylist_predictor );

	g_total_iteration += iterno;

	ylist[0] = yn;
	return yn;
}
