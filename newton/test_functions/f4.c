#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int nf = 1;

double x0[] = {
	1.5, 
};

void load_f ( double *x, double *f )
{
	double x1 = x[0];
	f[0] = x1*x1 - x1;
}


void load_jacobian ( double *x, double *J )
{
	double x1 = x[0];
	*(J + nf*0 + 0) = 2.0*x1 - 1.0; // (1,1)
}




