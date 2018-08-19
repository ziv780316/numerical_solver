#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "opts.h"

void show_help ();
void str_to_lower ( char *str );
int is_str_nocase_match ( const char *str_a, const char *str_b );

opt_t g_opts = {
	.iterative_type = NEWTON_NORMAL,
	.diff_type = NEWTON_DIFF_FORWARD,		
	.maxiter = -1, 	
	.rtol = 1e-3,	
	.atol = 1e-6,
	.residual_tol = 1e-9,
	.random_initial = false,
	.debug = false,
	.output_file = NULL
};

void show_help ()
{
	printf( "*------------------------------------*\n" 
		"*         Open Newton Solver         *\n"
		"*------------------------------------*\n" 
		"[Options]\n"
		"  -h  =>  show help\n"
		"  -d  =>  enable debug infomation\n"
		"  -z -->  randomize x0\n"
		"  -i | --iterative  =>  specify iterative method\n"
		"  -e | --derivative  =>  specify derivative type\n"
		"  -m | --maxiter  =>  specify maximum iterations\n"
		"  -r | --rtol  =>  specify rtol\n"
		"  -a | --atol  =>  specify atol\n"
		"  -u | --residual  =>  specify residual tol\n"
		"  -o | --output  =>  specify output file name\n"
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
			{"random_initial", no_argument, 0, 'z'},

			// setting options
			{"iterative", required_argument, 0, 'i'},
			{"derivative", required_argument, 0, 'e'},
			{"rtol", required_argument, 0, 'r'},
			{"atol", required_argument, 0, 'a'},
			{"residual", required_argument, 0, 'u'},
			{"output", required_argument, 0, 'o'},
			{"maxiter", required_argument, 0, 'm'},
			{0, 0, 0, 0}
		};

		// getopt_long stores the option index here
		int option_index = 0;

		c = getopt_long( argc, argv, "hdzi:e:r:a:t:o:m:u:", long_options, &option_index );

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

			case 'z':
				g_opts.random_initial = true;
				break;

			case 'o':
				g_opts.output_file = optarg;
				if ( !freopen( g_opts.output_file, "w", stdout ) )
				{
					fprintf( stderr, "[Error] open file fail --> %s\n", strerror(errno) );
					abort();
				}
				break;

			case 'i':
				if ( is_str_nocase_match( "normal", optarg ) )
				{
					g_opts.iterative_type = NEWTON_NORMAL;
				}
				else if ( is_str_nocase_match( "chord", optarg ) )
				{
					g_opts.iterative_type = NEWTON_CHORD;
				}
				else if ( is_str_nocase_match( "broyden", optarg ) )
				{
					g_opts.iterative_type = NEWTON_BROYDEN;
				}
				else if ( is_str_nocase_match( "broyden_inverted", optarg ) )
				{
					g_opts.iterative_type = NEWTON_BROYDEN_INVERTED;
				}
				else if ( is_str_nocase_match( "damped", optarg ) )
				{
					g_opts.iterative_type = NEWTON_DAMPED;
				}
				else if ( is_str_nocase_match( "line_search", optarg ) )
				{
					g_opts.iterative_type = NEWTON_LINE_SEARCH;
				}
				else
				{
					fprintf( stderr, "[Error] unknown newton iterative type %s\n", optarg );
					abort();
				}
				break;

			case 'e':
				if ( is_str_nocase_match( "jacobian", optarg ) )
				{
					g_opts.diff_type = NEWTON_DIFF_JACOBIAN;
				}
				else if ( is_str_nocase_match( "forward", optarg ) )
				{
					g_opts.diff_type = NEWTON_DIFF_FORWARD;
				}
				else if ( is_str_nocase_match( "central", optarg ) )
				{
					g_opts.diff_type = NEWTON_DIFF_CENTRAL;
				}
				else
				{
					fprintf( stderr, "[Error] unknown derivative approximation method %s\n", optarg );
					abort();
				}
				break;
				
			case 'r':
				g_opts.rtol = atof( optarg );
				break;

			case 'a':
				g_opts.atol = atof( optarg );
				break;

			case 'u':
				g_opts.residual_tol = atof( optarg );
				break;

			case 'm':
				g_opts.maxiter = atof( optarg );
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

