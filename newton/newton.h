#ifndef NEWTON_H
#define NEWTON_H

typedef enum {
	NEWTON_NORMAL,
	NEWTON_CHORD,
	NEWTON_JACOBI,
	NEWTON_BROYDEN,
	NEWTON_BROYDEN_INVERTED,
	NEWTON_BROYDEN_INVERTED_BAD,
} newton_iterative_type;

typedef enum {
	MODIFIED_NONE,
	MODIFIED_DAMPED,
	MODIFIED_LINE_SEARCH
} newton_modified_type;

typedef enum {
	RESCUE_NONE,
	RESCUE_DIAGONAL,
} newton_rescue_type;

typedef enum {
	NEWTON_DIFF_JACOBIAN,
	NEWTON_DIFF_FORWARD,
	NEWTON_DIFF_CENTRAL
} newton_derivative_type;

// perform Newton-Raphson iterations 
bool newton_solve ( newton_iterative_type iterative_type, 
		    newton_modified_type modified_type,
		    newton_rescue_type rescue_type,
		    newton_derivative_type diff_type,
		    int n,
		    double *x0, // initial x
		    double *x_result, // final x
		    double *f_result, // final f(x)
		    void (load_f) (double *x, double*f),
		    void (load_jacobian) (double *x, double*J),
		    int maxiter,
		    int miniter,
		    int *total_iter,
		    double rtol,
		    double atol,
		    double residual_tol,
		    double max_dx,
		    double jmin,
		    bool random_initial,
		    bool debug,
		    char *debug_file );

#endif

