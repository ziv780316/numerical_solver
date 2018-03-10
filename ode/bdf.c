#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "methods.h"

/*
definition:
	apply 'predictor-corrector method'

	GEAR-2
	yn_corrector = (4.0 / 3.0) * yn_1 - (1.0 / 3.0) * yn_2 + (2.0 / 3.0) * hn_1 * diff(yn_predictor, xn)
	hn_1 = xn - xn_1

reduce:	
	use 'explicit method' to predictor yn_predictor (forward-euler)

	yn_predictor = yn_1 + diff(yn_1, xn_1)
	yn_predictor = yn_1 + hn_1 * (yn_1 - yn_2) / hn_2 
*/
/*
double gear ( double xn_1, double hn_1, double *ylist, size_t order )
{
	double yn_1 = ylist[0];
	double yn_2 = ylist[1];
	double yn;
	double yn_predictor;
	double yn_corrector;
	double xn = xn_1 + hn_1;
	double tol = 1e-3;

	if ( 1 == order )
	{
		yn_predictor = forward_euler( xn_1, hn_1, ylist );
		yn_corrector = yn_1 + hn_1 * diff(yn_predictor, xn_1);

		while ( fabs(yn_corrector - yn_predictor) > tol )
		{
			yn_predictor = yn_corrector;
			yn_corrector = yn_1 + hn_1 * diff(yn_predictor, xn);
			//printf("%.10e %.10e\n", yn_corrector, yn_predictor );
		}
	}
	else if ( 2 == order )
	{
		yn_predictor = forward_euler( xn_1, hn_1, ylist );
		yn_corrector = (4.0 / 3.0) * yn_1 - (1.0 / 3.0) * yn_2 + (2.0 / 3.0) * hn_1 * diff(yn_predictor, xn);

		while ( fabs(yn_corrector - yn_predictor) > tol )
		{
			yn_predictor = yn_corrector;
			yn_corrector = (4.0 / 3.0) * yn_1 - (1.0 / 3.0) * yn_2 + (2.0 / 3.0) * hn_1 * diff(yn_predictor, xn);
			//printf("%.10e %.10e\n", yn_corrector, yn_predictor );
		}
	}
	else
	{
		double yn_3 = ylist[2];

		yn_predictor = forward_euler( xn_1, hn_1, ylist );
		yn_corrector = (18.0 / 11.0) * yn_1 - (9.0 / 11.0) * yn_2 + (2.0 / 11.0) * yn_3 + (6.0 / 11.0) * hn_1 * diff(yn_predictor, xn);

		while ( fabs(yn_corrector - yn_predictor) > tol )
		{
			yn_predictor = yn_corrector;
			yn_corrector = (18.0 / 11.0) * yn_1 - (9.0 / 11.0) * yn_2 + (2.0 / 11.0) * yn_3 + (6.0 / 11.0) * hn_1 * diff(yn_predictor, xn);
			//printf("%.10e %.10e\n", yn_corrector, yn_predictor );
		}
	}

	yn = yn_corrector;

	return yn;
}

*/
