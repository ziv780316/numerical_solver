#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int nf = 4;

double x0[] = {
	0.0, 
	0.0, 
	0.0, 
	0.0, 
};

// v1 1 0 1
// g1 1 2 1
// i2 2 0 0.1
// g2 2 3 1
// g3 3 0 1
void load_f ( double *x, double *f )
{	
	double g1 = 1;
	double g2 = 2;
	double g3 = 1;
	double i1 = 0.1;
	f[0] = x[1] - 1;
	f[1] = (x[1] - x[2]) * g1 + x[0];
	f[2] = (x[2] - x[1]) * g1 + (x[2] - x[3]) * g2 + i1;
	f[3] = (x[3] - x[2]) * g2 + x[3] * g3;
}


void load_jacobian ( double *x, double *J )
{
}

bool bypass_check ( double *x, double *f, double *dx ) 
{
	bool bypass_violate = true;
	return bypass_violate;
}

