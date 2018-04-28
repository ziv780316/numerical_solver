#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef IVP1
/*
stiff ODE problem:
	y' = -y
	y(0) = 1

exact solution:
	y = exp(-t)
*/
double g_initial_solution = 1.0;

double diff ( double y, double t )
{
	return -y;
}

double exact ( double t )
{
	return exp(-t);
}

double jacobian ( double y, double t )
{
	return -1.0;
}

#elif IVP2
/*
ODE:
	(vo - vi)/R = C*(dvo/dt)

	y' = -(y - vi)/(R*C)
	y(0) = 0

exact solution:
	y = 1 - exp(-t)
*/
double g_initial_solution = 0;

double diff ( double y, double t )
{
	double vi = 1.0;
	double r  = 1.0;
	double c  = 1.0;
	return -(y - vi)/(r*c);
}
double exact ( double t )
{
	return 1 - exp(-t);
}

#elif IVP3
/*
stiff ODE problem:
	y' = -5*(t+1)*y^2 + 5/(t+1) - 1/((t+1)^2)
	y(0) = 1

exact solution:
	y = 1/(t+1)
*/
double g_initial_solution = 1.0;

double diff ( double y, double t )
{
	return -5.0*(t+1.0)*(y*y) + 5.0/(t+1.0) - 1.0/((t+1.0)*(t+1.0));
}

double exact ( double t )
{
	return 1.0/(t+1.0);
}

#elif IVP4
/*
stiff ODE problem:
	y' = y*(1 - y)
	y(0) = 0.01

exact solution:
	y = e^t/(99+e^t)
*/
double g_initial_solution = 0.01;

double diff ( double y, double t )
{
	return y * (1 - y);
}

double exact ( double t )
{
	return exp(t) / (99.0 + exp(t));
}

double jacobian ( double y, double t )
{
	return 1 - (2.0 * y);
}

#endif
