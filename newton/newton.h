#ifndef NEWTON_H
#define NEWTON_H

typedef enum {
	NEWTON_NORMAL,
	NEWTON_CHORD,
	NEWTON_CHORD_WITH_BYPASS_CHECK,
	NEWTON_JACOBI,
	NEWTON_BROYDEN,
	NEWTON_BROYDEN_INVERTED,
	NEWTON_BROYDEN_INVERTED_BAD,
} newton_iterative_type;

typedef enum {
	DAMPED_NONE,
	DAMPED_DIRECT,
	DAMPED_LINE_SEARCH
} newton_damped_type;

typedef enum {
	RESCUE_NONE,
	RESCUE_DIAGONAL,
} newton_rescue_type;

typedef enum {
	NEWTON_DIFF_JACOBIAN,
	NEWTON_DIFF_FORWARD,
	NEWTON_DIFF_CENTRAL
} newton_derivative_type;

typedef struct 
{
	int n_iter;
	int n_mat_factor;
	int n_mat_solve;
	int n_f_load;
	int n_jac_load;
} performance_stat;

// perform Newton-Raphson iterations 
bool newton_solve ( newton_iterative_type iterative_type, 
		    newton_damped_type damped_type,
		    newton_rescue_type rescue_type,
		    newton_derivative_type diff_type,
		    int n,
		    double *x0, // initial x
		    double *x_result, // final x
		    double *f_result, // final f(x)
		    void (load_f) (double *x, double*f),
		    void (load_jacobian) (double *x, double*J),
		    bool (bypass_check) (double *x, double *f, double *dx),
		    int maxiter,
		    int miniter,
		    double rtol,
		    double atol,
		    double bypass_rtol,
		    double bypass_atol,
		    double residual_tol,
		    double max_dx,
		    double jmin,
		    bool random_initial,
		    performance_stat *nr_stat,
		    bool debug,
		    char *debug_file );

#endif

