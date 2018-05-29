#include <stdio.h>
#include <stdlib.h>
#include "methods.h"

double rk ( int order, double tn_1, double yn_1, double hn_1, double *ylist )
{
	double yn;
	double f1 = hn_1 * diff(yn_1, tn_1);
	double f2 = hn_1 * diff(yn_1 + 0.5*f1, tn_1 + 0.5*hn_1);
	double f3 = hn_1 * diff(yn_1 + 0.5*f2, tn_1 + 0.5*hn_1);;
	double f4 = hn_1 * diff(yn_1 + f3, tn_1 + hn_1);;

	switch ( order )
	{
		case 4:
			yn = yn_1 + (1.0/6.0) * (f1 + 2.0*f2 + 2.0*f3 + f4);
			break;
		default:
			fprintf( stderr, "[Error] cannot support RK order=%d\n", order );
			abort();
			break;
	}

	return yn;
}

