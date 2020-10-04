#ifndef OPTS_H
#define OPTS_H

#include <stdbool.h>

#include "homotopy.h"

typedef struct
{
	homotopy_param_t homotopy_param;
	newton_param_t newton_param;
	char *output_file;
	char *problem_so;
	char *initial_x0_file;
} opt_t;

extern void show_help ();
extern void parse_cmd_options ( int argc, char **argv );
extern opt_t g_opts;

// user defined function
extern int nf;
extern double lamda;
extern double x0[];
extern void load_f ( double *x, double *y );
extern void load_jacobian ( double *x, double *J );
extern void bypass_check ( double *x, double *f, double *dx );

#endif
