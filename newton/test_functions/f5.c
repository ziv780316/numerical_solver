#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int nf = 2;

double x0[] = {
	2,
	1
};

void load_f ( double *x, double *f )
{
	double x1 = x[0];
	double x2 = x[1];
	f[0] = x1*x1 - x2*x2;
	f[1] = 3.0*x1*x1 - 3.0*x2*x2;
}


void load_jacobian ( double *x, double *J )
{
	double x1 = x[0];
	double x2 = x[1];
	*(J + nf*0 + 0) = 2.0*x1; // (1,1)
	*(J + nf*0 + 1) = 6.0*x1; // (2,1)
	*(J + nf*1 + 0) = -2.0*x2; // (1,2)
	*(J + nf*1 + 1) = -6.0*x2; // (2,2)
}




