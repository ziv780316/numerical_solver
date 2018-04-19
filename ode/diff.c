#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef IVP1
/*
stiff ODE problem:
	y' = -y
	y(0) = 1

exact solution:
	y = exp(-x)
*/
double g_initial_solution = 1.0;

double diff ( double y, double x )
{
	return -y;
}

double exact ( double x )
{
	return exp(-x);
}

#elif IVP2
/*
ODE:
	(vo - vi)/R = C*(dvo/dt)

	y' = -(y - vi)/(R*C)
	y(0) = 0

exact solution:
	y = 1 - exp(-x)
*/
double g_initial_solution = 0;

double diff ( double y, double x )
{
	double vi = 1.0;
	double r  = 1.0;
	double c  = 1.0;
	return -(y - vi)/(r*c);
}
double exact ( double x )
{
	return 1 - exp(-x);
}

#elif IVP3
/*
stiff ODE problem:
	y' = -5*(x+1)*y^2 + 5/(x+1) - 1/((x+1)^2)
	y(0) = 1

exact solution:
	y = 1/(x+1)
*/
double g_initial_solution = 1.0;

double diff ( double y, double x )
{
	return -5.0*(x+1.0)*(y*y) + 5.0/(x+1.0) - 1.0/((x+1.0)*(x+1.0));
}

double exact ( double x )
{
	return 1.0/(x+1.0);
}
#endif
