#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "matrix_solver.h"

void load_f ( double *x, double *f );
void load_jacobian ( double *x, double *J );
void load_df_dp ( double *x, double *df_dp );

int nf = 1;

double x0[] = {
	-1,
};

double p0 = 0;
double p = 0;

// p(x) = -x² + 1
// f(x,p) = x² - 1 + p = 0
void load_f ( double *x, double *f )
{
	double x1 = x[0];
	f[0] = x1*x1 - 1 + p;
}


void load_jacobian ( double *x, double *J )
{
	double x1 = x[0];
	*(J + nf*0 + 0) = 2.0 * x1; // (1,1)
}

void load_df_dp ( double *x, double *df_dp )
{
	double x1 = x[0];
	df_dp[0] = 1.0;
}

