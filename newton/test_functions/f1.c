#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int nf = 2;

double x0[] = {
	2.0, 
	2.0
};

void load_f ( double *x, double *f )
{
	f[0] = x[0] + cos(x[1]);
	f[1] = x[1] + sin(x[0]);
}


void load_jacobian ( double *x, double *J )
{
	*(J + nf*0 + 0) = 1.0; // (1,1)
	*(J + nf*0 + 1) = cos(x[0]); // (2,1)
	*(J + nf*1 + 0) = -sin(x[1]); // (1,2)
	*(J + nf*1 + 1) = 1.0; // (2,2)
}


