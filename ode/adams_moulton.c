#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "methods.h"

/*
definition:
	apply 'predictor-corrector method'

	yn_corrector = yn_1 + hn_1 * diff(yn_predict, xn)
	hn_1 = xn - xn_1

reduce:	
	use 'explicit method' to predictor yn_predictor (forward-euler)

	yn_predictor = yn_1 + diff(yn_1, xn_1)
	yn_predictor = yn_1 + hn_1 * (yn_1 - yn_2) / hn_2 
*/
/*
double backward_euler ( double xn_1, double hn_1, double *ylist )
{
	double yn_1 = ylist[0];
	double yn;
	double yn_diff;
	double xn = xn_1 + hn_1;
	double tol = 1e-3;

	// predictor-corrector method (explicit method, not ensure implicit method stability)
	double yn_predictor;
	double yn_corrector;
	yn_predictor = forward_euler( xn_1, hn_1, ylist );
	//yn_corrector = yn_1 + hn_1 * diff(yn_predictor, xn_1);
	//printf("xx %.10e %.10e\n", yn_corrector, yn_predictor );

	//while ( fabs(yn_corrector - yn_predictor) > tol )
	//{
	//	yn_predictor = yn_corrector;
	//	yn_corrector = yn_1 + hn_1 * diff(yn_predictor, xn);
	//	//printf("%.10e %.10e\n", yn_corrector, yn_predictor );
	//}

	//yn = yn_corrector;

	// Newton-Raphoson method
	//yn = yn_predictor; // NR initial value of this time point
	yn = yn_1; // NR initial value of this time point
	int iterno = 0;
	double nr_diff;
	double yn_next;
	double g;
	double g_delta;
	double g_diff;
	double delta = 1e-6;
	do 
	{
		++iterno;

		// use forward difference approximation
		g = yn - yn_1 - hn_1 * diff(yn, xn);
		g_delta = (yn + delta) - yn_1 - hn_1 * diff((yn + delta), xn);
		g_diff = (g_delta - g) / delta;
		yn_next = yn - (g / g_diff);
		nr_diff = fabs(yn_next - yn);
		//printf( "%d: yn_next=%.10e yn=%.10e yn_1=%.10e g=%.10e, g_delta=%.10e, g_diff=%.10e\n", iterno, yn_next, yn, yn_1, g, g_delta, g_diff );
		yn = yn_next;
	} while ( nr_diff > tol );
	

	return yn;
}
*/

double adams_moulton ( int order, double xn_1, double yn_1, double hn_1, double *ylist )
{
	double yn;
	double diff_1 = diff(yn_1, xn_1);
	double diff_2 = ylist[1];
	double diff_3 = ylist[2];
	double diff_4 = ylist[3];
	ylist[0] = diff_1;

	switch ( order )
	{
		case 1:
			// Backward-Euler
			yn = yn_1 + hn_1 * diff_1;
			break;
		case 2:
			yn = yn_1 + hn_1 * (
			  (3.0/2.0) * diff_1
			- (1.0/2.0) * diff_2 
			);
			break;
		case 3:
			yn = yn_1 + hn_1 * (
			  (23.0/12.0) * diff_1
			- (4.0/3.0) * diff_2 
			+ (5.0/12.0) * diff_3
			);
			break;
		case 4:
			yn = yn_1 + hn_1 * (
			  (55.0/24.0) * diff_1
			- (59.0/24.0) * diff_2 
			+ (37.0/24.0) * diff_3
			- (3.0/8.0) * diff_4 
			);
			break;
		default:
			fprintf( stderr, "[Error] cannot support AB order=%d\n", order );
			abort();
			break;
	}

	return yn;
}
