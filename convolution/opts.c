#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include "opts.h"

opt_t g_opts = {
	.freq_response_file = NULL,
	.input_signal_time_domain_file = NULL,
	.output_signal_time_domain_file = NULL,
	.configure_file = NULL,
};

void show_help ()
{
	printf( "*------------------------------------*\n" 
		"*         Convoluton Solver          *\n"
		"*------------------------------------*\n" 
		"[Convoluton Options]\n"
		"  -h  =>  show help\n"
		"  -d  =>  debug\n"
		"  -t | --time  =>  specify input signal time domain file\n"
		"  -f | --freq  =>  specify freq response file\n"
		"  -o | --output  =>  specify convolution result file\n"
		"  -c | --configure  =>  specify configure file\n"
		);
}

void parse_cmd_options ( int argc, char **argv )
{
	int c;

	while ( true )
	{
		static struct option long_options[] =
		{
			// flag options
			{"help", no_argument, 0, 'h'},
			{"debug", no_argument, 0, 'd'},

			// setting options
			{"time", required_argument, 0, 't'},
			{"freq", required_argument, 0, 'f'},
			{"output", required_argument, 0, 'o'},
			{"configure", required_argument, 0, 'c'},

			{0, 0, 0, 0}
		};

		// getopt_long stores the option index here
		int option_index = 0;

		c = getopt_long( argc, argv, "hdt:f:o:c:", long_options, &option_index );

		// detect the end of the options
		if ( -1 == c )
		{
			break;
		}

		switch ( c )
		{
			case 'h':
				show_help();
				exit( EXIT_SUCCESS );
				break;

			case 't':
				g_opts.input_signal_time_domain_file = optarg;
				break;

			case 'f':
				g_opts.freq_response_file = optarg;
				break;

			case 'o':
				g_opts.output_signal_time_domain_file = optarg;
				break;

			case 'c':
				g_opts.configure_file = optarg;
				break;

			case 'd':
				g_opts.debug = true;
				break;

			case '?':
				/* getopt_long already printed an error message. */
				break;

			default:
				abort ();
				break;
		}
	}

	// print any remaining command line arguments (not options)
	if (optind < argc)
	{
		fprintf( stderr, "[Warning] non-option ARGV-elements: " );
		while ( optind < argc )
		{
			fprintf( stderr, "%s ", argv[optind++] );
		}
		fprintf( stderr, "\n" );
	}
}

