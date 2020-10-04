#ifndef HOMOTOPY_H
#define HOMOTOPY_H

#include "newton.h"

typedef enum {
	HOMOTOPY_EXTRAPOLATE_NONE,
	HOMOTOPY_EXTRAPOLATE_DIFFERENCE,
	HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL,
} homotopy_extrapolate_type;

typedef enum {
	HOMOTOPY_DF_DP_EXACT,
	HOMOTOPY_DF_DP_FORWARD,
	HOMOTOPY_DF_DP_CENTRAL
} homotopy_df_dp_type;

typedef struct 
{
	int n_step;
	int n_success;
	int n_fail;
	int n_limit_point;

	int n_iter;
	int n_mat_factor;
	int n_mat_solve;
	int n_mat_solve_sensitivity;
	int n_f_load;
	int n_f_load_sensitivity;
	int n_df_dp_load;
	int n_jac_load;
} homotopy_performance_stat_t;

typedef struct
{
	homotopy_extrapolate_type extrapolate_type;
	homotopy_df_dp_type df_dp_type;

	double lamda_start;
	double lamda_stop;
	double arc_length;

	bool debug;

	homotopy_performance_stat_t hom_stat;

} homotopy_param_t;

// perform BBD Newton-Raphson iterations 
bool arc_length_bbd_newton_solve ( 
		newton_param_t *newton_param,
		int *perm, // permuation for matrix ordering
		double *J, // jacobian
		double *pp,
		double *x0, // initial x
	        double *x_ans, // user give solution x*
		double p0, // initial p
		double *x_result, // final x
		double *p_result, // final x
		double *f_result, // final f(x)
		double *g_result, // final g(x)
		double *xc, // hyper circle central (xc,pc)
		double pc, // hyper circle central (xc,pc)
		double dt, // arc length
		void (load_f) (double *x, double*f),
		void (load_df_dp) (double *x, double*f),
		void (load_jacobian) (double *x, double*J),
		char *debug_file );


#endif

