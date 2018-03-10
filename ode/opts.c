#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include "opts.h"

void show_help ();
void str_to_lower ( char *str );
int is_str_nocase_match ( const char *str_a, const char *str_b );

opt_t g_opts = {
	.method = BDF,	
	.solver = NEWTON,		
	.tstep = 1e-3,
	.tstop = 100,
	.iteraton_limit = -1, // unlimit
	.maxord = 1,
	.debug = false
};

void show_help ()
{
	printf( "*------------------------------------*\n" 
		"*           Open ODE Solver          *\n"
		"*------------------------------------*\n" 
		"[Options]\n"
		"  -h  =>  show help\n"
		"  -d  =>  enable debug infomation\n"
		"  -t | --tstep [num]  =>  specify integration size\n"
		"  -p | --tstop [num]  =>  specify stop time\n"
		"  -m | --method \"[name]\"  =>  specify integration method\n"
		"  -s | --solver \"[name]\"  =>  specify solver type\n"
		"  -o | --order \"[num]\"  =>  specify integration order\n"
		"  -i | --iteration \"[num]\"  =>  specify iteration limit each time\n"
		);
}

void str_to_lower ( char *str )
{
	for ( int i = 0; '\0' != str[i]; ++i )
	{
		str[i] = tolower( str[i] );
	}
}

int is_str_nocase_match ( const char *str_a, const char *str_b )
{
	char *a = (char *) calloc ( strlen(str_a) + 1, sizeof(char) );
	char *b = (char *) calloc ( strlen(str_b) + 1, sizeof(char) );
	strcpy( a, str_a );
	strcpy( b, str_b );
	str_to_lower( a );
	str_to_lower( b );
	return (0 == strcmp( a, b ));
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
			{"tstep", required_argument, 0, 't'},
			{"tstop", required_argument, 0, 'p'},
			{"method", required_argument, 0, 'm'},
			{"solver", required_argument, 0, 's'},
			{"order", required_argument, 0, 'o'},
			{"iteration_limit", required_argument, 0, 'i'},
			{0, 0, 0, 0}
		};

		// getopt_long stores the option index here
		int option_index = 0;

		c = getopt_long( argc, argv, "hdm:s:o:t:p:i:", long_options, &option_index );

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

			case 'd':
				g_opts.debug = true;
				break;

			case 't':
				g_opts.tstep = atof( optarg );
				break;

			case 'p':
				g_opts.tstop = atof( optarg );
				break;
				
			case 's':
				if ( is_str_nocase_match( optarg, "newton" ) )
				{
					g_opts.method = NEWTON;
				}
				else if ( is_str_nocase_match( optarg, "predictor_corrector" ) )
				{
					g_opts.method = PREDICTOR_CORRECTOR;
				}
				else
				{
					fprintf( stderr, "[Error] unknow integration method '%s'\n", optarg );
					abort();
				}
				break;

			case 'o':
				g_opts.maxord = atoi( optarg );
				break;

			case 'i':
				g_opts.iteraton_limit = atoi( optarg );
				break;

			case 'm':
				if ( is_str_nocase_match( optarg, "am" ) )
				{
					g_opts.method = AM;
				}
				else if ( is_str_nocase_match( optarg, "ab" ) )
				{
					g_opts.method = AB;
				}
				else if ( is_str_nocase_match( optarg, "bdf" ) )
				{
					g_opts.method = BDF;
				}
				else if ( is_str_nocase_match( optarg, "simpson" ) )
				{
					g_opts.method = SIMPSON;
				}
				else if ( is_str_nocase_match( optarg, "rk" ) )
				{
					g_opts.method = RK;
				}
				else
				{
					fprintf( stderr, "[Error] unknow integration method '%s'\n", optarg );
					abort();
				}
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
