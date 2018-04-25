#include <stdio.h>
#include <stdlib.h>
#include "methods.h"

double adams_bashforth ( int order, double tn_1, double yn_1, double hn_1, double *ylist )
{
	double yn;
	double diff_1 = diff(yn_1, tn_1);
	double diff_2 = ylist[1];
	double diff_3 = ylist[2];
	double diff_4 = ylist[3];
	ylist[0] = diff_1;

	switch ( order )
	{
		case 1:
			// Forward-Euler
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

