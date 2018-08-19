#ifndef NEWTON_H
#define NEWTON_H

typedef enum {
	NEWTON_NORMAL,
	NEWTON_CHORD,
	NEWTON_BROYDEN,
	NEWTON_BROYDEN_INVERTED,
	NEWTON_DAMPED,
	NEWTON_LINE_SEARCH
} newton_iterative_type;

typedef enum {
	NEWTON_DIFF_JACOBIAN,
	NEWTON_DIFF_FORWARD,
	NEWTON_DIFF_CENTRAL
} newton_derivative_type;

#define NEWTON_DIFF_FORWARD 0
#define NEWTON_DIFF_CENTRAL 1
#define NEWTON_DIFF_JACOBIAN 2

// perform newton iterations 
bool newton_solve ( newton_iterative_type iterative_type, newton_derivative_type diff_type, int maxiter, double rtol, double atol, double residual_tol, bool random_initial, bool debug );

// define by package user
void load_f ( double *x, double *y ) __attribute__((weak));
void load_jacobian ( double *x, double *J ) __attribute__((weak)); // J is column-major matrix
extern int n_f;
extern double x0[];

#endif

