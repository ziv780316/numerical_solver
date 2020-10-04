#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "matrix_solver.h"

void load_f ( double *x, double *J );
void load_jacobian ( double *x, double *J );
void load_df_dp ( double *x, double *df_dp );

int nf = 1;

double x0[] = {
	-1,
};

double p = 0;

void load_f ( double *x, double *f )
{
	double x1 = x[0];
	double c = (1 - p);
	f[0] = (x1 - c) * (x1 - c);
}


void load_jacobian ( double *x, double *J )
{
	double x1 = x[0];
	double c = (1 - p);
	*(J + nf*0 + 0) = 2.0 * (x1 - c); // (1,1)
}

void load_df_dp ( double *x, double *df_dp )
{
	double x1 = x[0];
	double c = (1 - p);
	df_dp[0] = -2.0 * (x1 - c);
}
