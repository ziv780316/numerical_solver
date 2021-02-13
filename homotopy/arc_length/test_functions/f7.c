#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "matrix_solver.h"

void load_f ( double *x, double *f );
void load_jacobian ( double *x, double *J );
void load_df_dp ( double *x, double *df_dp );

int nf = 2;

double x0[] = {
	0.5,
	0.5,
};

double p0 = 0;
double p = 0;

// p(x) = -1.5*x³ - 1.5*x² + 0.5*x + 0.8
// f1(x1,x2,p) = 1.5*x1³ + 0.2*x2³+ 1.5*x1² - 0.5*x1 - 0.8 + p = 0
// f2(x1,x2,p) = 1.4*x2³ + 0.1*x1³ + 1.6*x2² - 0.4*x2 - 0.7 + p = 0
void load_f ( double *x, double *f )
{
	double x1 = x[0];
	double x2 = x[1];
	f[0] = 1.5*x1*x1*x1 + 0.2*x2*x2*x2 + 1.5*x1*x1 - 0.5*x1 - 0.8 + p;
	f[1] = 1.4*x2*x2*x2 + 0.1*x1*x1*x1 + 1.6*x2*x2 - 0.4*x2 - 0.7 + p;
}


void load_jacobian ( double *x, double *J )
{
	double x1 = x[0];
	double x2 = x[1];
	*(J + nf*0 + 0) = 4.5*x1*x1 + 3*x1 - 0.5; // (1,1)
	*(J + nf*0 + 1) = 0.6*x2*x2; // (1,2)
	*(J + nf*1 + 0) = 0.3*x1*x1; // (2,1)
	*(J + nf*1 + 1) = 4.2*x2*x2 + 3.2*x2 - 0.4; // (2,2)
}

void load_df_dp ( double *x, double *df_dp )
{
	double x1 = x[0];
	double x2 = x[1];
	df_dp[0] = 1.0;
	df_dp[1] = 1.0;
}

