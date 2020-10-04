#ifndef HOMOTOPY_H
#define HOMOTOPY_H

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

	bool debug;

	homotopy_performance_stat_t hom_stat;

} homotopy_param_t;


#endif

