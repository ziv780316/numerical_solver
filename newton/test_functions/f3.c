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

void load_f ( double *x, double *f )
{
	double x1 = x[0];
	double x2 = x[1];
	double x3  = x[2];
	f[0] = 1 + 10*sin(x1) + sin(x2) + sin(x3);
	f[1] = 2 + cos(x1) + 10*cos(x2) + cos(x3);
	f[2] = 2 + x1*x2 + 10*x3*x3;
}


void load_jacobian ( double *x, double *J )
{
	double x1 = x[0];
	double x2 = x[1];
	double x3  = x[2];
	*(J + nf*0 + 0) = 10*cos(x1); // (1,1)
	*(J + nf*0 + 1) = -sin(x1); // (2,1)
	*(J + nf*0 + 2) = x2; // (3,1)
	*(J + nf*1 + 0) = cos(x2); // (1,2)
	*(J + nf*1 + 1) = -10*sin(x2); // (2,2)
	*(J + nf*1 + 2) = x1; // (3,2)
	*(J + nf*2 + 0) = cos(x3); // (1,3)
	*(J + nf*2 + 1) = -sin(x3); // (2,3)
	*(J + nf*2 + 2) = 20*x3; // (3,3)
}



