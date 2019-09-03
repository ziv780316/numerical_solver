#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include "opts.h"

static void str_to_lower ( char *str );
static int is_str_nocase_match ( const char *str_a, const char *str_b );

opt_t g_opts = {
	.iterative_type = NEWTON_NORMAL,
	.damped_type = DAMPED_NONE,
	.rescue_type = RESCUE_NONE,
	.diff_type = NEWTON_DIFF_FORWARD,		
	.maxiter = -1, 	
	.miniter = 0, 	
	.rtol = 1e-3,	
	.atol = 1e-6,
	.bypass_rtol = 1e-2,	
	.bypass_atol = 1e-3,
	.max_dx = DBL_MAX,
	.jmin = 0.0,
	.residual_tol = 1e-9,
	.random_initial = false,
	.debug = false,
	.output_file = NULL,
	.problem_so = NULL,
	.initial_x0_file = NULL
};

void show_help ()
{
	printf( "*------------------------------------*\n" 
		"*         Open Newton Solver         *\n"
		"*------------------------------------*\n" 
		"[Options]\n"
		"  -h  =>  show help\n"
		"  -d  =>  enable debug information\n"
		"  -z  =>  randomize x0\n"
		"  -i | --iterative  =>  specify iterative method (default normal)\n"
		"    normal\n"
		"    chord\n"
		"    chord_bypass_check\n"
		"    jacobi\n"
		"    broyden\n"
		"    broyden_inverted\n"
		"    broyden_inverted_bad\n"
		"  -f | --damped  =>  specify damped type (default none)\n"
		"    damped\n"
		"    line_search\n"
		"  -c | --rescue  =>  specify rescue method (default none)\n"
		"    diagonal\n"
		"  -e | --derivative  =>  specify derivative type (default forward)\n"
		"    jacobian\n"
		"    forward\n"
		"    central\n"
		"  -m | --maxiter  =>  specify maximum iterations (default unlimited)\n"
		"  -n | --miniter  =>  specify minimum iterations (default unlimited)\n"
		"  -r | --rtol  =>  specify rtol (default 1e-3)\n"
		"  -a | --atol  =>  specify atol (default 1e-6)\n"
		"  -y | --bypass rtol  =>  specify bypass_rtol (default 1e-2)\n"
		"  -s | --bypass atol  =>  specify bypass_atol (default 1e-3)\n"
		"  -u | --residual  =>  specify residual tol (default 1e-9)\n"
		"  -b | --max_dx  =>  maximum dx in damped newton (default unlimited)\n"
		"  -j | --jmin  =>  specify minimum diagonal value of jacobian (default 0)\n"
		"  -o | --output  =>  specify output file name (default terminal)\n"
		"  -p | --problem_so  =>  specify problem file (*.so)\n"
		"  -x | --initial_x0_file  =>  specify initializing file of x0\n"
		);
}

static void str_to_lower ( char *str )
{
	for ( int i = 0; '\0' != str[i]; ++i )
	{
		str[i] = tolower( str[i] );
	}
}

static int is_str_nocase_match ( const char *str_a, const char *str_b )
{
	char *a = (char *) calloc ( strlen(str_a) + 1, sizeof(char) );
	char *b = (char *) calloc ( strlen(str_b) + 1, sizeof(char) );
	bool is_same;
	strcpy( a, str_a );
	strcpy( b, str_b );
	str_to_lower( a );
	str_to_lower( b );
	is_same = (0 == strcmp( a, b ));
	free( a );
	free( b );
	return is_same;
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
			{"damped", required_argument, 0, 'f'},
			{"rescue", required_argument, 0, 'c'},
			{"derivative", required_argument, 0, 'e'},
			{"rtol", required_argument, 0, 'r'},
			{"atol", required_argument, 0, 'a'},
			{"bypass_rtol", required_argument, 0, 'y'},
			{"bypass_atol", required_argument, 0, 's'},
			{"residual", required_argument, 0, 'u'},
			{"max_dx", required_argument, 0, 'b'},
			{"jmin", required_argument, 0, 'j'},
			{"output", required_argument, 0, 'o'},
			{"maxiter", required_argument, 0, 'm'},
			{"miniter", required_argument, 0, 'n'},
			{"problem_so", required_argument, 0, 'p'},
			{"initial_x0_file", required_argument, 0, 'x'},
			{0, 0, 0, 0}
		};

		// getopt_long stores the option index here
		int option_index = 0;

		c = getopt_long( argc, argv, "hdzi:f:c:e:r:a:y:s:t:o:m:n:u:b:j:p:x:", long_options, &option_index );

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
				else if ( is_str_nocase_match( "jacobi", optarg ) )
				{
					g_opts.iterative_type = NEWTON_JACOBI;
				}
				else if ( is_str_nocase_match( "chord", optarg ) )
				{
					g_opts.iterative_type = NEWTON_CHORD;
				}
				else if ( is_str_nocase_match( "chord_bypass_check", optarg ) )
				{
					g_opts.iterative_type = NEWTON_CHORD_WITH_BYPASS_CHECK;
				}
				else if ( is_str_nocase_match( "broyden", optarg ) )
				{
					g_opts.iterative_type = NEWTON_BROYDEN;
				}
				else if ( is_str_nocase_match( "broyden_inverted", optarg ) )
				{
					g_opts.iterative_type = NEWTON_BROYDEN_INVERTED;
				}
				else if ( is_str_nocase_match( "broyden_inverted_bad", optarg ) )
				{
					g_opts.iterative_type = NEWTON_BROYDEN_INVERTED_BAD;
				}
				else
				{
					fprintf( stderr, "[Error] unknown newton iterative type %s\n", optarg );
					abort();
				}
				break;

			case 'f':
				if ( is_str_nocase_match( "damped", optarg ) )
				{
					g_opts.damped_type = DAMPED_DIRECT;
				}
				else if ( is_str_nocase_match( "line_search", optarg ) )
				{
					g_opts.damped_type = DAMPED_LINE_SEARCH;
				}
				else
				{
					fprintf( stderr, "[Error] unknown damped type %s\n", optarg );
					abort();
				}
				break;

			case 'c':
				if ( is_str_nocase_match( "diagonal", optarg ) )
				{
					g_opts.rescue_type = RESCUE_DIAGONAL;
				}
				else
				{
					fprintf( stderr, "[Error] unknown rescue type %s\n", optarg );
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

			case 'y':
				g_opts.bypass_rtol = atof( optarg );
				break;

			case 's':
				g_opts.bypass_atol = atof( optarg );
				break;

			case 'u':
				g_opts.residual_tol = atof( optarg );
				break;

			case 'b':
				g_opts.max_dx = atof( optarg );
				break;

			case 'j':
				g_opts.jmin = atof( optarg );
				break;

			case 'm':
				g_opts.maxiter = atoi( optarg );
				break;

			case 'n':
				g_opts.miniter = atoi( optarg );
				break;

			case 'p':
				g_opts.problem_so = optarg;
				break;

			case 'x':
				g_opts.initial_x0_file = optarg;

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

