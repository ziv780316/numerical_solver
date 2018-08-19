#ifndef OPTS_H
#define OPTS_H

#include <stdbool.h>

#include "newton.h"

typedef struct
{
	newton_iterative_type iterative_type;		
	newton_derivative_type diff_type;		
	int maxiter;	
	double rtol;	
	double atol;	
	double residual_tol;	
	bool random_initial;
	bool debug;
	char *output_file;
} opt_t;

extern void parse_cmd_options ( int argc, char **argv );
extern opt_t g_opts;

#endif
