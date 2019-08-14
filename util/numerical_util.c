#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "numerical_util.h"

void FUNC_WITH_POSTFIX( interpolation_lagrange, POSTFIX ) (
	long n, // Pn(xp) = yp
	long ldy, // length of yp
	double *x, 
	double **y,
	double xp,
	double *yp
)
{
	if ( n < 0 )
	{
		fprintf( stderr, "[Error] interpolation order %ld should >= 0\n", n );
		abort();
	}

	FLOAT_TYPE *l = (FLOAT_TYPE *) malloc ( sizeof(FLOAT_TYPE) * (n + 1) );

	for ( long i = 0; i <= n; ++i )
	{
		l[i] = 1.0;
		for ( long k = 0; k <= n; ++k )
		{
			if ( i != k )
			{
				l[i] *= (xp - x[k]) / (x[i] - x[k]);
			}
		}
	}

	for ( long i = 0; i < ldy; ++i )
	{
		yp[i] = 0.0;
		for ( long k = 0; k <= n; ++k )
		{
			yp[i] += y[i][k] * l[k];
		}
	}

	free( l );
}



