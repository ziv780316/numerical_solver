#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <strings.h>

#include "numerical_util.h"

void FUNC_WITH_POSTFIX( interpolation_lagrange, POSTFIX ) (
	long n, // Pn(xp) = yp
	long ldy, // length of yp
	double *x, 
	double **y, // backward y0, y1, y2 --> y[n], y[n-1]...,y[0]
	double xp,
	double *yp
)
{
	if ( n < 0 )
	{
		fprintf( stderr, "[Error] interpolation order %ld should >= 0\n", n );
		abort();
	}
	else if ( 0 == n )
	{
		bcopy( yp, y[0], sizeof(double) * ldy );
		return;
	}

	FLOAT_TYPE *l = (FLOAT_TYPE *) malloc ( sizeof(FLOAT_TYPE) * (n + 1) );

	FLOAT_TYPE _xp = xp;
	FLOAT_TYPE _xk;
	FLOAT_TYPE _xi;
	for ( long i = 0; i <= n; ++i )
	{
		l[i] = 1.0;
		_xi = x[i];
		for ( long k = 0; k <= n; ++k )
		{
			_xk = x[k];
			if ( i != k )
			{
				l[i] *= (_xp - _xk) / (_xi - _xk);
			}
		}
	}

	FLOAT_TYPE _yik;
	FLOAT_TYPE _ypi;
	for ( long i = 0; i < ldy; ++i )
	{
		_ypi = 0.0;
		for ( long k = 0; k <= n; ++k )
		{
			_yik = y[i][k];
			_ypi += _yik * l[k];
		}
		yp[i] = _ypi;
	}

	free( l );
}


// backward difference y0 = y[0], yn = y[n]
// ddN = dd[n][0], dd1 = dd[1][0]
// ddN = [yn, yn-1, ..., y0] = (y⁽ⁿ⁾(ζ) / N!) where xn ≤ ζ ≤ x0
void FUNC_WITH_POSTFIX( divide_difference, POSTFIX ) (
	long n, // backward difference [yn, yn-1, ..., y0], i.e. (y1 - y0)/(x1 - x0) and x0 > x1
	double *x, // x0 is x[0], xn is x[n]
	double *y, // backward difference, y0 is y[0], yn is y[n]
	double *dd 
)
{
	if ( n <= 0 )
	{
		fprintf( stderr, "[Error] divide difference order %ld should > 0\n", n );
		abort();
	}

	FLOAT_TYPE *_dd = (FLOAT_TYPE *) malloc ( sizeof(FLOAT_TYPE) * (n + 1) );
	for ( long i = 0; i <= n; ++i )
	{
		_dd[i] = y[i];
	}


	for ( long i = 1; i <= n; ++i )
	{
		FLOAT_TYPE _xk;
		FLOAT_TYPE _xki;
		for ( long k = n; k >= i; --k )
		{
			// bottom up approach
			// i = 1
			// dd[1] = [y1, y0]
			//       = (y[1] - y[0]) / (x[1] - x[0]) 
			//       = (dd[1] - dd[0]) / (x[1] - x[0]) 
			//
			// i = 2
			// dd[2] = [y2, y1, y0]
			//       = [y2, y1] - [y1, y0] 
			//       = (dd[2] - dd[1]) / (x[2] - x[0]) 
			//
			// i = n
			// dd[n] = [yn, yn-1, ..., y0] 
			//       = [yn, yn-1, ..., 1] - [yn-1, yn-2, ..., 0]
			//       = (dd[n] - dd[n-1]) / (x[n] - x[0]) 
			_xk = x[k];
			_xki = x[k - i];
			_dd[k] = (_dd[k] - _dd[k - 1]) / (_xk - _xki); 
		}
	}

	for ( long i = 0; i <= n; ++i )
	{
		dd[i] = _dd[i];
	}
	free( _dd );
}

