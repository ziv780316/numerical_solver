#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "opts.h"
#include "methods.h"

double adams_moulton ( int order, double xn_1, double yn_1, double hn_1, double *ylist )
{
	double xn = xn_1 + hn_1;
	double yn;
	double yn_predictor;
	double yn_k;
	double yn_k_delta;
	double diff_0;
	double diff_0_delta;
	double diff_1 = ylist[0];
	double diff_2 = ylist[1];
	double diff_3 = ylist[2];

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
		yn_predictor = adams_bashforth( order, xn_1, yn_1, hn_1, ylist_predictor );
		yn = yn_predictor; // NR initial value with predictor
	}
	else
	{
		yn = yn_1; // NR initial value without predictor
	}

	do 
	{
		diff_0       = diff(yn, xn);
		diff_0_delta = diff(yn + delta, xn); // forward difference

		switch ( order )
		{
			case 1:
				// Backward-Euler
				yn_k = yn_1 + hn_1 * diff_0;
				yn_k_delta = yn_1 + hn_1 * diff_0_delta;
				break;
			case 2:
				// Trapezoidal
				yn_k = yn_1 + hn_1 * (1.0/2.0) * (diff_0 + diff_1);
				yn_k_delta = yn_1 + hn_1 * (1.0/2.0) * (diff_0_delta + diff_1);
				break;
			case 3:
				yn_k = yn_1 + hn_1 * (
					  (5.0/12.0) * diff_0
					+ (2.0/3.0) * diff_1
					- (1.0/12.0) * diff_2
					);
				yn_k_delta = yn_1 + hn_1 * (
					  (5.0/12.0) * diff_0_delta
					+ (2.0/3.0) * diff_1
					- (1.0/12.0) * diff_2
					);
				break;
			case 4:
				yn_k = yn_1 + hn_1 * (
					  (3.0/8.0) * diff_0
					+ (19.0/24.0) * diff_1
					- (5.0/24.0) * diff_2
					+ (1.0/24.0) * diff_3
					);
				yn_k_delta = yn_1 + hn_1 * (
					  (3.0/8.0) * diff_0_delta
					+ (19.0/24.0) * diff_1
					- (5.0/24.0) * diff_2
					+ (1.0/24.0) * diff_3
					);
				break;
			default:
				fprintf( stderr, "[Error] cannot support AM order=%d\n", order );
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

	ylist[0] = diff_0;
	free( ylist_predictor );

	g_total_iteration += iterno;

	return yn;
}
