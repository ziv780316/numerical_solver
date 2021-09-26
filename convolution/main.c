#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "opts.h"
#include "conv.h"
#include "c_python.h"

int main ( int argc, char **argv )
{
	if ( 1 == argc )
	{
		show_help();
	}
	else
	{
		// getopt parse command line arguments
		parse_cmd_options ( argc, argv );

		// read frequency response
		if ( !g_opts.freq_response_file )
		{
			fprintf( stderr, "[Error] must specify freq response file by -f\n" );
			exit(1);
		}
		FILE *fin = fopen( g_opts.freq_response_file, "r" );
		if ( !fin )
		{
			fprintf( stderr, "[Error] fopen '%s' fail --> %s\n", g_opts.freq_response_file, strerror(errno) );
			exit(1);
		}

		conv_db_t *conv_db = (conv_db_t *) calloc ( 1, sizeof(conv_db_t) );
		int lineno = 0;
		double f;
		double re;
		double im;
		char buf[BUFSIZ];
		while ( fgets(buf, BUFSIZ, fin ) )
		{
			++lineno;
			int n_read = sscanf( buf, "%le %le %le\n", &f, &re, &im );
			if ( n_read != 3 )
			{
				fprintf( stderr, "[Warning] line=%d read fail\n + string='%s'\n", lineno, buf );
			}
			else
			{
				++(conv_db->n_freq_point);
				conv_db->f = (double *) realloc ( conv_db->f, sizeof(double)*conv_db->n_freq_point );
				conv_db->H = (complex *) realloc ( conv_db->H, sizeof(complex)*conv_db->n_freq_point );
				(conv_db->f)[conv_db->n_freq_point - 1] = f;
				(conv_db->H)[conv_db->n_freq_point - 1] = re + I*im;
			}
		}
		if ( g_opts.debug )
		{
			printf( "[DEBUG] fstart=%.10le fend=%.10le n_freq_point=%d\n", (conv_db->f)[0], (conv_db->f[conv_db->n_freq_point-1]), conv_db->n_freq_point );
		}

		// read time domain signal 
		if ( !g_opts.input_signal_time_domain_file )
		{
			fprintf( stderr, "[Error] must specify input time domain file by -t\n" );
			exit(1);
		}
		fin = fopen( g_opts.input_signal_time_domain_file, "r" );
		if ( !fin )
		{
			fprintf( stderr, "[Error] fopen '%s' fail --> %s\n", g_opts.input_signal_time_domain_file, strerror(errno) );
			exit(1);
		}

		lineno = 0;
		double t;
		double x;
		while ( fgets(buf, BUFSIZ, fin ) )
		{
			++lineno;
			int n_read = sscanf( buf, "%le %le\n", &t, &x );
			if ( n_read != 2 )
			{
				fprintf( stderr, "[Warning] line=%d read fail\n + string='%s'\n", lineno, buf );
			}
			else
			{
				++(conv_db->n_time_point);
				conv_db->t = (double *) realloc ( conv_db->t, sizeof(double)*conv_db->n_time_point );
				conv_db->x = (double *) realloc ( conv_db->x, sizeof(double)*conv_db->n_time_point );
				(conv_db->t)[conv_db->n_time_point - 1] = t;
				(conv_db->x)[conv_db->n_time_point - 1] = x;
			}
		}
		if ( g_opts.debug )
		{
			printf( "[DEBUG] tstart=%.10le tend=%.10le n_time_point=%d\n", (conv_db->t)[0], (conv_db->t[conv_db->n_time_point-1]), conv_db->n_time_point );
		}

		// read configure
		if ( !g_opts.configure_file )
		{
			fprintf( stderr, "[Warning] there is no configure file specified, use default convolution settings\n" );
			conv_init_default_setting( conv_db );
		}
		else
		{
			fin = fopen( g_opts.configure_file, "r" );
			if ( !fin )
			{
				fprintf( stderr, "[Error] fopen '%s' fail --> %s\n", g_opts.configure_file, strerror(errno) );
				exit(1);
			}
		}
		conv_db->debug = g_opts.debug;
		conv_evaluate_basic( conv_db );
		conv_show_configure( conv_db );

		// evaluate h by performing IFFT on H 
		conv_get_impulse_response ( conv_db );
	}
	
	return EXIT_SUCCESS;
}

