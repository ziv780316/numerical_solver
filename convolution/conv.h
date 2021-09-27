#ifndef CONV_H
#define CONV_H

#include <complex.h>

typedef struct
{
	// raw data
	int n_freq_point;
	double *f;
	double *f_interp;
	complex *H;
	complex *H_interp;
	double *H_mag_interp;

	int n_time_point;
	double *t;
	double *x;

	// impulse response by ifft
	double *t_interp;
	double *h;

	// convolution settings
	double fmax;
	double fdelta;
	double fs;
	double ts;
	double T;

	int H_dc_extrap_method;
	#define CONV_H_DC_EXTRAP_CONSTANT 0
	#define CONV_H_DC_EXTRAP_LINEAR 1

	int H_hf_extrap_method;
	#define CONV_H_HF_EXTRAP_CONSTANT 0
	#define CONV_H_HF_EXTRAP_PHASE_LINEAR 1

	int n_sample; // one-side
	int n; // 2 * n_sample
	int H_interp_method;
	#define CONV_H_INTERP_LINEAR 0
	#define CONV_H_INTERP_SPLINE 1

	double h_truncate_tail_ratio;
	int h_padding_type;
	#define CONV_H_ZERO_PADDING 0
	#define CONV_H_IDENTICAL_REPEAT_PADDING 1
	#define CONV_H_DECAY_REPEAT_PADDING 2

	int integration_method;
	#define CONV_INTEGRATION_BE 0
	#define CONV_INTEGRATION_TRAP 1

	bool debug;

} conv_db_t;

void conv_init_default_setting ( conv_db_t *conv_db );
void conv_evaluate_basic ( conv_db_t *conv_db );
void conv_show_configure ( conv_db_t *conv_db );
void conv_interp_H ( conv_db_t *conv_db );
void conv_get_impulse_response ( conv_db_t *conv_db );

#endif

