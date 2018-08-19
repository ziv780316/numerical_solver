#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int n_f = 2;

double x0[] = {2.0, 2.0};
//double x0[] = {-0.9, 1.0};

void load_f ( double *x, double *y )
{
	y[0] = x[0] + cos(x[1]);
	y[1] = x[1] + sin(x[0]);
}


void load_jacobian ( double *x, double *J )
{
	*(J + n_f*0 + 0) = 1.0; // (1,1)
	*(J + n_f*0 + 1) = cos(x[0]); // (2,1)
	*(J + n_f*1 + 0) = -sin(x[1]); // (1,2)
	*(J + n_f*1 + 1) = 1.0; // (2,2)
}


