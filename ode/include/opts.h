#ifndef OPTS_H
#define OPTS_H

#include <stdbool.h>

typedef enum
{
	AM, // Adams-Moulton
	AB, // Adams-Bashforth
	BDF, // Backward Differential Formula
	SIMPSON, // Simpson's method
	RK // Runge-Kutta
} INTEGRATION_TYPE;

typedef enum
{
	PREDICTOR_CORRECTOR,
	NEWTON
} SOLVER_TYPE;

typedef struct
{
	INTEGRATION_TYPE method;		
	SOLVER_TYPE solver;		
	double tstep;	
	double tstop;	
	int iteraton_limit;
	int maxord;
	bool debug;
} opt_t;

extern void parse_cmd_options ( int argc, char **argv );

#endif
