#ifndef OPTS_H
#define OPTS_H

#include <stdbool.h>

#include "newton.h"

typedef struct
{
	newton_iterative_type iterative_type;		
	newton_damped_type damped_type;		
	newton_rescue_type rescue_type;		
	newton_derivative_type diff_type;		
	int maxiter;	
	int miniter;	
	double rtol;	
	double atol;	
	double bypass_rtol;	
	double bypass_atol;	
	double residual_tol;	
	double max_dx;
	double jmin;
	bool random_initial;
	bool debug;
	char *output_file;
	char *problem_so;
	char *initial_x0_file;
} opt_t;

extern void show_help ();
extern void parse_cmd_options ( int argc, char **argv );
extern opt_t g_opts;

// user defined function
extern int nf;
extern double x0[];
extern void load_f ( double *x, double *y );
extern void load_jacobian ( double *x, double *J );
extern void bypass_check ( double *x, double *f, double *dx );

#endif
