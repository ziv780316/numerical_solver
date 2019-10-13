#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int nf = 3;

double x0[] = {
	1.0, 
	1.0,
	1.0
};

double x_ans[] = {
	7.071067811865476e-01,
	7.071067811865476e-01,
	7.071067811865476e-01
};


void load_f ( double *x, double *f )
{
	double x1 = x[0];
	double x2 = x[1];
	double v  = x[2];
	f[0] = 1 - (2.0 * v * x1);
	f[1] = 1 - (2.0 * v * x2);
	f[2] = (x1 * x1) + (x2 * x2) - 1.0;
}


void load_jacobian ( double *x, double *J )
{
	double x1 = x[0];
	double x2 = x[1];
	double v  = x[2];
	*(J + nf*0 + 0) = -2.0 * v; // (1,1)
	*(J + nf*0 + 1) = 0.0; // (2,1)
	*(J + nf*0 + 2) = 2.0 * x1; // (3,1)
	*(J + nf*1 + 0) = 0.0; // (1,2)
	*(J + nf*1 + 1) = -2.0 * v; // (2,2)
	*(J + nf*1 + 2) = 2.0 * x2; // (3,2)
	*(J + nf*2 + 0) = -2.0 * x1; // (1,3)
	*(J + nf*2 + 1) = -2.0 * x2; // (2,3)
	*(J + nf*2 + 2) = 0.0; // (3,3)
}



