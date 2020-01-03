#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

int nf = 3;

double x0[] = {
	0.0, 
	0.0, 
	0.0, 
};

// v1 1 0 1
//
// Y matrix
// ---------
// g1 1 2 1 
// g2 2 3 1 
// ---------
// 
// g3 3 0 1
void load_f ( double *x, double *f )
{	
	double g1 = 1;
	double g2 = 1;
	double g3 = 1;
	double y11 = g1*g2/(g1 + g2);
	double y21 = -y11;
	double y12 = -y11;
	double y22 = y11;
	f[0] = x[1] - 1;
	f[1] = x[1] * y11 + x[2] * y12 + x[0];
	f[2] = x[1] * y21 + x[2] * y22 + x[2] * g3;
}


void load_jacobian ( double *x, double *J )
{
}

bool bypass_check ( double *x, double *f, double *dx ) 
{
	bool bypass_violate = true;
	return bypass_violate;
}

