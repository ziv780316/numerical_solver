#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <complex.h>

// v = dc
void generate_v_dc ( char *output, double tstart, double tend, double tdelta, double dc );

// H = Vo/Vin = 1/(1+1j*w*r*c)
void generate_freq_rc_lpf ( char *output, double fstart, double fend, double fdelta, double r, double c );

int main ( int argc, char **argv )
{
	generate_v_dc( "vdc_time_signal", 0, 10, 0.1, 1 );
	generate_freq_rc_lpf( "rc_freq_response", 0, 1e3, 0.1, 1, 1 );
	
	return EXIT_SUCCESS;
}

// --------------------------------------------------------------------------
// H = Vo/Vin = 1/(1+1j*w*r*c)
void generate_freq_rc_lpf ( char *output, double fstart, double fend, double fdelta, double r, double c )
{
	FILE *fout = fopen( output, "w" );
	if ( !fout )
	{
		fprintf( stderr, "[Error] fopen '%s' fail --> %s\n", output, strerror(errno) );
		exit(1);
	}

	long double f = fstart;
	complex H;
	for ( int i = 0; ; ++i )
	{
		f = fdelta * i;

		H = 1 / (1 + I*2*M_PI*f*r*c);
		fprintf( fout, "%.15Le %.15le %.15le\n", f, creal(H), cimag(H) );
		if ( ((fdelta * (i+1)) > fend) && ((fend - f) > 1e-3 * fdelta) )
		{
			f = fend;
		}

		if ( f >= fend )
		{
			break;
		}
	}
	fclose( fout );
}

void generate_v_dc ( char *output, double tstart, double tend, double tdelta, double dc )
{
	FILE *fout = fopen( output, "w" );
	if ( !fout )
	{
		fprintf( stderr, "[Error] fopen '%s' fail --> %s\n", output, strerror(errno) );
		exit(1);
	}

	long double t = 0;
	double v;
	for ( int i = 0; ; ++i )
	{
		t = tdelta * i;
		v = dc;
		fprintf( fout, "%.15Le %.15le\n", t, v );
		if ( ((tdelta * (i+1)) > tend) && ((tend - t) > 1e-3 * tdelta) )
		{
			t = tend;
		}

		if ( t >= tend )
		{
			break;
		}
	}
	fclose( fout );
}
