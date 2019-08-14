#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "numerical_util.h"

double f1 ( double x ) 
{
	return x;
}

double f2 ( double x ) 
{
	return x * x;
}

double f3 ( double x ) 
{
	return x * x * x; 
}

int main ( int argc, char **argv )
{
	long n = 2;
	long ldy = 3;
	double x[3] = {3, 1, 0};
	double **y;
	double xp = 2;
	double yp[3];
	y = (double **) malloc ( ldy * sizeof(double *) );
	for ( int i = 0; i < ldy; ++i )
	{
		y[i] = (double *) malloc ( n * sizeof(double) );
	}
	y[0][0] = f1(x[0]);
	y[0][1] = f1(x[1]);
	y[0][2] = f1(x[2]);
	y[1][0] = f2(x[0]);
	y[1][1] = f2(x[1]);
	y[1][2] = f2(x[2]);
	y[2][0] = f3(x[0]);
	y[2][1] = f3(x[1]);
	y[2][2] = f3(x[2]);

	FUNC_WITH_POSTFIX( interpolation_lagrange, POSTFIX ) ( 1, ldy, x, y, xp, yp );
	for ( long i = 0; i < ldy; ++i )
	{
		for ( long k = 0; k <= n; ++k )
		{
			printf( "y%ld[%ld] = %.15le\n", i, k, y[i][k] );
		}
		printf( "yp[%ld] = %.15le\n\n", i, yp[i] );
	}

	return EXIT_SUCCESS;
}

