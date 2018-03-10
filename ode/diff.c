#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
stiff ODE problem:
	y' = -y
	y(0) = 1

exact solution:
	y = exp(-x)
*/
double g_initial_solution = 1.0;

double diff ( double yn, double xn )
{
	if ( 0.0 == yn )
	{
		return -1.0;
	}
	return -yn;
}

double diff2 ( double yn, double xn )
{
	return -1.0;
}

double diff_exact ( double yn, double xn )
{
	return -exp( -xn );
}


double exact ( double xn )
{
	return exp( -xn );
}

/*
ODE:
	(vo - vi) / R = C * (dvo / dt)

	y' = -(y - vi) / (R * C)
	y(0) = 0

exact solution:
	y = 1 - exp(-x)
*/
//double g_initial_solution = 0;
//
//double diff ( double yn, double xn )
//{
//	double vi = 1.0;
//	double r  = 1.0;
//	double c  = 1.0;
//	return -(yn - vi) / (r * c);
//}
//double exact ( double xn )
//{
//	return 1 - exp( -xn );
//}
