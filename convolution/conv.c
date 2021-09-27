#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "conv.h"
#include "fftw3.h"
#include "c_python.h"

void conv_init_default_setting ( conv_db_t *conv_db )
{
	conv_db->H_hf_extrap_method = CONV_H_HF_EXTRAP_PHASE_LINEAR;
	conv_db->H_dc_extrap_method = CONV_H_DC_EXTRAP_CONSTANT;
	conv_db->n_sample = 4096;
	conv_db->h_truncate_tail_ratio = 1.0;
	conv_db->h_padding_type = CONV_H_ZERO_PADDING;
	conv_db->integration_method = CONV_INTEGRATION_BE;
	conv_db->debug = false;
}

void conv_evaluate_basic ( conv_db_t *conv_db )
{
	if ( 0 == (conv_db->n_sample % 2) )
	{
		conv_db->n= 2 * conv_db->n_sample;
		conv_db->fmax = (conv_db->f)[conv_db->n_freq_point - 1];
		conv_db->fdelta = conv_db->fmax / conv_db->n_sample;
		conv_db->T = 1 / conv_db->fdelta;
		conv_db->fs = 2 * conv_db->fmax; 
		conv_db->ts = 1 / conv_db->fs;
	}
	else
	{
		fprintf( stderr, "[Error] only support even n_sample\n" );
		exit(1);
	}
}

void conv_show_configure ( conv_db_t *conv_db )
{
	char *str;
	printf( "----------------------------------------------\n" );
	printf( "* Convoluton configure:\n" );
	printf( "\n" );
	printf( " + fmax = %.15le\n", conv_db->fmax );
	printf( " + n_sample = %d\n", conv_db->n_sample );
	printf( " + fdelta = %.15le\n", conv_db->fdelta );
	printf( " + fs = %.15le\n", conv_db->fs );
	printf( " + T = %.15le\n", conv_db->T );
	printf( " + ts = %.15le\n", conv_db->ts );
	printf( " + h_truncate_tail_ratio = %.15le\n", conv_db->h_truncate_tail_ratio );
	switch ( conv_db->H_dc_extrap_method )
	{
		case CONV_H_DC_EXTRAP_CONSTANT: str = "constant"; break;
		case CONV_H_DC_EXTRAP_LINEAR: str = "linear"; break;
		default: fprintf( stderr, "[Error] unknown H_dc_extrap_method=%d\n", conv_db->H_dc_extrap_method ); exit(1); break;
	}
	printf( " + H_dc_interp_method = %s\n", str );
	switch ( conv_db->H_hf_extrap_method )
	{
		case CONV_H_HF_EXTRAP_CONSTANT: str = "constant"; break;
		case CONV_H_HF_EXTRAP_PHASE_LINEAR: str = "phase_linear"; break;
		default: fprintf( stderr, "[Error] unknown H_hf_extrap_method=%d\n", conv_db->H_hf_extrap_method ); exit(1); break;
	}
	printf( " + H_hf_interp_method = %s\n", str );
	switch ( conv_db->H_interp_method )
	{
		case CONV_H_INTERP_LINEAR: str = "linear"; break;
		case CONV_H_INTERP_SPLINE: str = "spline"; break;
		default: fprintf( stderr, "[Error] unknown H_interp_method=%d\n", conv_db->H_interp_method ); exit(1); break;
	}
	printf( " + H_interp_method = %s\n", str );
	switch ( conv_db->h_padding_type )
	{
		case CONV_H_ZERO_PADDING: str = "zero"; break;
		case CONV_H_IDENTICAL_REPEAT_PADDING: str = "identical"; break;
		case CONV_H_DECAY_REPEAT_PADDING: str = "decay"; break;
		default: fprintf( stderr, "[Error] unknown h_padding_type=%d\n", conv_db->h_padding_type ); exit(1); break;
	}
	printf( " + h_padding_type = %s\n", str );
	switch ( conv_db->integration_method )
	{
		case CONV_INTEGRATION_BE: str = "euler"; break;
		case CONV_INTEGRATION_TRAP: str = "trape"; break;
		case CONV_H_DECAY_REPEAT_PADDING: str = "decay_repeat"; break;
		default: fprintf( stderr, "[Error] unknown integration_method=%d\n", conv_db->integration_method ); exit(1); break;
	}
	printf( " + integration_method = %s\n", str );
	printf( "\n" );
	printf( "----------------------------------------------\n" );
}

void conv_interp_H (conv_db_t *conv_db )
{
	if ( CONV_H_INTERP_LINEAR == conv_db->H_interp_method )
	{
		conv_db->f_interp = (double *) calloc ( conv_db->n, sizeof(double) );
		conv_db->H_interp = (complex *) calloc ( conv_db->n, sizeof(complex) );
		conv_db->H_mag_interp = (double *) calloc ( conv_db->n, sizeof(double) );
		bool find_interval;
		int k = 1;
		complex Hp;
		double f_left;
		double f_right;
		double mag;
		double mag_left;
		double mag_right;
		double phase;
		double phase_left;
		double phase_right;
		double fp;
		double fstart = (conv_db->f)[0];
		double fend = (conv_db->f)[conv_db->n_freq_point - 1];
		double slope;
		double df;
		double fdelta = conv_db->fmax / conv_db->n_sample;
		double f_interval;
		for ( int i = 0; i <= conv_db->n_sample; ++i )
		{
			fp = fdelta * i;
			if ( i == 0 )
			{
				// DC extrapolation
				if ( 0 == fstart )
				{
					Hp = (conv_db->H)[0];
				}
				else
				{
					if ( CONV_H_DC_EXTRAP_CONSTANT == conv_db->H_dc_extrap_method )
					{
						mag = cabs( (conv_db->H)[0] );
						Hp = mag;
						if ( creal((conv_db->H)[0]) < 0 )
						{
							Hp *= -1;
						}
					}
					else
					{
						fprintf( stderr, "[Error] unknown H_dc_extrap_method=%d\n", conv_db->H_dc_extrap_method );
						exit(1);
					}
				}
			}
			else if ( fp > fend )
			{
				// HF extrapolation
				if ( CONV_H_HF_EXTRAP_PHASE_LINEAR == conv_db->H_hf_extrap_method )
				{
					f_left = conv_db->f[conv_db->n_freq_point - 2];
					f_right = conv_db->f[conv_db->n_freq_point - 1];
					mag = cabs( (conv_db->H)[conv_db->n_freq_point - 1] );
					phase_left = carg( conv_db->H[conv_db->n_freq_point - 2] );
					phase_right = carg( conv_db->H[conv_db->n_freq_point - 1] );
					f_interval = f_right - f_left;
					slope = (phase_right - phase_left) / f_interval;
					df = fp - f_left;
					phase = phase_left + (df * slope);
					Hp = (mag * cos(phase)) + I*(mag * sin(phase));
				}
				else
				{
					fprintf( stderr, "[Error] unknown H_dc_extrap_method=%d\n", conv_db->H_dc_extrap_method );
					exit(1);
				}
			}
			else
			{
				// interpolation
				if ( fp < fstart )
				{
					find_interval = true;
					f_left = 0;
					f_right = conv_db->f[0];
					mag_left = cabs((conv_db->H_interp)[0]);
					if ( creal((conv_db->H_interp)[0]) >= 0 )
					{
						phase_left = 0;
					}
					else
					{
						phase_left = M_PI;
					}
					mag_right = cabs((conv_db->H)[0]);
					phase_right = carg((conv_db->H)[0]);
				}
				else
				{
					find_interval = false;
					for ( ; k < conv_db->n_freq_point; ++k )
					{
						if ( (fp <= (conv_db->f)[k]) && (fp > (conv_db->f)[k - 1]) )
						{
							find_interval = true;
							f_left = conv_db->f[k - 1];
							f_right = conv_db->f[k];
							mag_left = cabs((conv_db->H)[k - 1]);
							phase_left = carg((conv_db->H)[k - 1]);
							mag_right = cabs((conv_db->H)[k]);
							phase_right = carg((conv_db->H)[k]);
							break;
						}
					}
				}

				if ( !find_interval )
				{
					fprintf( stderr, "[Error] assert (find_interval == true) in %s lineno=%d\n", __FILE__, __LINE__ );
					fprintf( stderr, " + fp=%.10le k=%d\n", fp, k );
					exit(1);
				}

				df = fp - f_left;
				f_interval = f_right - f_left;
				slope = (mag_right - mag_left) / f_interval;
				mag = mag_left + (df * slope);
				slope = (phase_right - phase_left) / f_interval;
				phase = phase_left + (df * slope);
				Hp = (mag * cos(phase)) + I*(mag * sin(phase));
			}

			conv_db->f_interp[i] = fp;
			conv_db->H_interp[i] = Hp;
			conv_db->H_mag_interp[i] = cabs(Hp);
		}
	}
	else
	{
		fprintf( stderr, "[Error] unknown H_interp_method=%d\n", conv_db->H_interp_method );
		exit(1);
	}

	if ( conv_db->debug )
	{
		printf( "----------------------------------------------\n" );
		printf( "* H_interp:\n" );
		printf( "\n" );
		for ( int i = 0; i <= conv_db->n_sample; ++i )
		{
			printf( "%d: %.15le %.15le %.15le\n", i+1, (conv_db->f_interp)[i], creal((conv_db->H_interp)[i]), cimag((conv_db->H_interp)[i]) );
		}
		printf( "\n" );
		printf( "----------------------------------------------\n" );
	}
}

void conv_get_impulse_response ( conv_db_t *conv_db )
{
	// interp
	conv_interp_H( conv_db );

	//// ifft
	unsigned flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT;
	conv_db->t_interp = (double *) malloc ( sizeof(double) * conv_db->n);
	conv_db->h = (double *) malloc ( sizeof(double) * conv_db->n );
	fftw_plan plan = fftw_plan_dft_c2r_1d( conv_db->n,  // 2 * n_sample
	                                       conv_db->H_interp, // only store positive H_interp of length (n/2)+1 from 0 to fmax
					       conv_db->h, 
					       flags ); 
	fftw_execute( plan );
	fftw_destroy_plan( plan );

	for ( int i = 0; i < conv_db->n; ++i )
	{
		(conv_db->t_interp)[i] = conv_db->ts * i;
		(conv_db->h)[i] *= conv_db->fdelta;
	}

	if ( conv_db->debug )
	{
		printf( "----------------------------------------------\n" );
		printf( "* Impulse response of IFFT H_interp:\n" );
		printf( "\n" );
		for ( int i = 0; i < conv_db->n; ++i )
		{
			printf( "%.15le %.15le\n", (conv_db->t_interp)[i], (conv_db->h)[i] );
		}
		printf( "\n" );
		printf( "----------------------------------------------\n" );

		python_init( 0 ); // debug = 0
		python_create_list( "t", conv_db->t_interp, conv_db->n, C_PYTHON_VALUE_TYPE_FLOAT_64 );
		python_create_list( "h", conv_db->h, conv_db->n, C_PYTHON_VALUE_TYPE_FLOAT_64 );
		python_create_list( "f", conv_db->f_interp, conv_db->n_sample + 1, C_PYTHON_VALUE_TYPE_FLOAT_64 );
		python_create_list( "H", conv_db->H_mag_interp, conv_db->n_sample + 1, C_PYTHON_VALUE_TYPE_FLOAT_64 );
		python_eval_string( 
			"fig, axs = plt.subplots(2,2,figsize=(12,8));" // figsize unit is monitor inches
			"fig.suptitle('impulse response h(t)');"
			"axs[0,0].plot(f,H,'bo-',markersize=3);"
			"axs[0,0].set_ylabel('H(f)');"
			"axs[0,0].set_xlabel('Hz');"
			"axs[0,0].set_xscale('log');"
			"axs[0,1].plot(t,h,'b-',markersize=1);"
			"axs[0,1].set_ylabel('h(t)');"
			"axs[0,1].set_xlabel('t');"
			"plt.show();" 
		);
		python_close();
	}
}
