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
	.homotopy_param = 
	{
		.extrapolate_type = HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL,
		.df_dp_type = HOMOTOPY_DF_DP_FORWARD,
		.lamda_start = 0,
		.lamda_stop = 1,
		.debug = false,
		.hom_stat = {0}
	},
	.newton_param = 
	{
		.iterative_type = NEWTON_NORMAL,
		.damped_type = DAMPED_NONE,
		.rescue_type = RESCUE_NONE,
		.diff_type = NEWTON_DIFF_FORWARD,		
		.n = -1,
		.maxiter = -1, 	
		.miniter = 0, 	
		.delta_rtol = 1e-3,	
		.delta_atol = 1e-6,
		.bypass_rtol = 1e-2,	
		.bypass_atol = 1e-3,
		.max_dx = DBL_MAX,
		.jmin = 0.0,
		.line_search_tol = 1e-3,
		.residual_rtol = 1e-3,
		.residual_atol = 1e-9,
		.random_initial = false,
		.debug = false,
		.nr_stat = {0}
	},
	.output_file = NULL,
	.problem_so = NULL,
	.initial_x0_file = NULL
};

void show_help ()
{
	printf( "*------------------------------------*\n" 
		"*       Open Homotopy Solver         *\n"
		"*------------------------------------*\n" 
		"[Newton Options]\n"
		"  -h  =>  show help\n"
		"  -z  =>  randomize x0\n"
		"  -d | --debug =>  enable debug information (default none)\n"
		"    newton\n"
		"    homotopy\n"
		"    all\n"
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
		"  -e | --derivative  =>  specify ∂f/∂x type (default forward)\n"
		"    jacobian\n"
		"    forward\n"
		"    central\n"
		"  -m | --maxiter  =>  specify maximum iterations (default unlimited)\n"
		"  -n | --miniter  =>  specify minimum iterations (default unlimited)\n"
		"  -r | --delta_rtol  =>  specify delta_rtol (default 1e-3)\n"
		"  -a | --delta_atol  =>  specify delta_atol (default 1e-6)\n"
		"  -g | --residual_rtol  =>  specify residual rtol (default 1e-3)\n"
		"  -u | --residual_atol  =>  specify residual atol (default 1e-9)\n"
		"  -y | --bypass rtol  =>  specify bypass_rtol (default 1e-2)\n"
		"  -s | --bypass atol  =>  specify bypass_atol (default 1e-3)\n"
		"  -b | --max_dx  =>  maximum dx in damped newton (default unlimited)\n"
		"  -j | --jmin  =>  specify minimum diagonal value of jacobian (default 0)\n"
		"  -o | --output  =>  specify output file name (default terminal)\n"
		"  -p | --problem_so  =>  specify problem file (*.so)\n"
		"  -x | --initial_x0_file  =>  specify initializing file of x0\n"
		"\n"
		"[Homotopy Options]\n"
		"  -t | --extrapolation_type  =>  specify extrapolate type (default none)\n"
		"    difference\n"
		"    differential\n"
		"  -k | --derivative  =>  specify ∂f/∂p type (default forward)\n"
		"    exact\n"
		"    forward\n"
		"    central\n"
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
			{"random_initial", no_argument, 0, 'z'},

			// setting options
			{"debug", required_argument, 0, 'd'},
			{"iterative", required_argument, 0, 'i'},
			{"damped", required_argument, 0, 'f'},
			{"rescue", required_argument, 0, 'c'},
			{"derivative", required_argument, 0, 'e'},
			{"delta_rtol", required_argument, 0, 'r'},
			{"delta_atol", required_argument, 0, 'a'},
			{"residual_rtol", required_argument, 0, 'g'},
			{"residual_atol", required_argument, 0, 'u'},
			{"bypass_rtol", required_argument, 0, 'y'},
			{"bypass_atol", required_argument, 0, 's'},
			{"max_dx", required_argument, 0, 'b'},
			{"jmin", required_argument, 0, 'j'},
			{"line_search_tol", required_argument, 0, 'l'},
			{"output", required_argument, 0, 'o'},
			{"maxiter", required_argument, 0, 'm'},
			{"miniter", required_argument, 0, 'n'},
			{"problem_so", required_argument, 0, 'p'},
			{"initial_x0_file", required_argument, 0, 'x'},

			// homotopy options
			{"extrapolate_type", required_argument, 0, 't'},
			{"df_dp_derivative", required_argument, 0, 'k'},
			{0, 0, 0, 0}
		};

		// getopt_long stores the option index here
		int option_index = 0;

		c = getopt_long( argc, argv, "hzd:i:f:c:e:r:a:y:s:t:o:m:n:g:u:b:j:l:p:x:t:k:", long_options, &option_index );

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
				if ( is_str_nocase_match( "newton", optarg ) )
				{
					g_opts.newton_param.debug = true;
				}
				else if ( is_str_nocase_match( "homotopy", optarg ) )
				{
					g_opts.homotopy_param.debug = true;
				}
				else if ( is_str_nocase_match( "all", optarg ) )
				{
					g_opts.newton_param.debug = true;
					g_opts.homotopy_param.debug = true;
				}
				break;

			case 'z':
				g_opts.newton_param.random_initial = true;
				break;

			case 'i':
				if ( is_str_nocase_match( "normal", optarg ) )
				{
					g_opts.newton_param.iterative_type = NEWTON_NORMAL;
				}
				else if ( is_str_nocase_match( "jacobi", optarg ) )
				{
					g_opts.newton_param.iterative_type = NEWTON_JACOBI;
				}
				else if ( is_str_nocase_match( "chord", optarg ) )
				{
					g_opts.newton_param.iterative_type = NEWTON_CHORD;
				}
				else if ( is_str_nocase_match( "chord_bypass_check", optarg ) )
				{
					g_opts.newton_param.iterative_type = NEWTON_CHORD_WITH_BYPASS_CHECK;
				}
				else if ( is_str_nocase_match( "broyden", optarg ) )
				{
					g_opts.newton_param.iterative_type = NEWTON_BROYDEN;
				}
				else if ( is_str_nocase_match( "broyden_inverted", optarg ) )
				{
					g_opts.newton_param.iterative_type = NEWTON_BROYDEN_INVERTED;
				}
				else if ( is_str_nocase_match( "broyden_inverted_bad", optarg ) )
				{
					g_opts.newton_param.iterative_type = NEWTON_BROYDEN_INVERTED_BAD;
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
					g_opts.newton_param.damped_type = DAMPED_DIRECT;
				}
				else if ( is_str_nocase_match( "line_search", optarg ) )
				{
					g_opts.newton_param.damped_type = DAMPED_LINE_SEARCH;
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
					g_opts.newton_param.rescue_type = RESCUE_DIAGONAL;
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
					g_opts.newton_param.diff_type = NEWTON_DIFF_JACOBIAN;
				}
				else if ( is_str_nocase_match( "forward", optarg ) )
				{
					g_opts.newton_param.diff_type = NEWTON_DIFF_FORWARD;
				}
				else if ( is_str_nocase_match( "central", optarg ) )
				{
					g_opts.newton_param.diff_type = NEWTON_DIFF_CENTRAL;
				}
				else
				{
					fprintf( stderr, "[Error] unknown derivative approximation method %s\n", optarg );
					abort();
				}
				break;
				
			case 'r':
				g_opts.newton_param.delta_rtol = atof( optarg );
				break;

			case 'a':
				g_opts.newton_param.delta_atol = atof( optarg );
				break;

			case 'y':
				g_opts.newton_param.bypass_rtol = atof( optarg );
				break;

			case 's':
				g_opts.newton_param.bypass_atol = atof( optarg );
				break;

			case 'g':
				g_opts.newton_param.residual_rtol = atof( optarg );
				break;

			case 'u':
				g_opts.newton_param.residual_atol = atof( optarg );
				break;

			case 'b':
				g_opts.newton_param.max_dx = atof( optarg );
				break;

			case 'j':
				g_opts.newton_param.jmin = atof( optarg );
				break;

			case 'l':
				g_opts.newton_param.line_search_tol = atof( optarg );
				break;

			case 'm':
				g_opts.newton_param.maxiter = atoi( optarg );
				break;

			case 'n':
				g_opts.newton_param.miniter = atoi( optarg );
				break;

			case 'o':
				g_opts.output_file = optarg;
				if ( !freopen( g_opts.output_file, "w", stdout ) )
				{
					fprintf( stderr, "[Error] open file fail --> %s\n", strerror(errno) );
					abort();
				}
				break;

			case 'p':
				g_opts.problem_so = optarg;
				break;

			case 'x':
				g_opts.initial_x0_file = optarg;
				break;

			case 't':
				if ( is_str_nocase_match( "none", optarg ) )
				{
					g_opts.homotopy_param.extrapolate_type = HOMOTOPY_EXTRAPOLATE_NONE;
				}
				else if ( is_str_nocase_match( "difference", optarg ) )
				{
					g_opts.homotopy_param.extrapolate_type = HOMOTOPY_EXTRAPOLATE_DIFFERENCE;
				}
				else if ( is_str_nocase_match( "differential", optarg ) )
				{
					g_opts.homotopy_param.extrapolate_type = HOMOTOPY_EXTRAPOLATE_DIFFERENTIAL;
				}
				break;

			case 'k':
				if ( is_str_nocase_match( "exact", optarg ) )
				{
					g_opts.homotopy_param.df_dp_type = HOMOTOPY_DF_DP_EXACT;
				} else if ( is_str_nocase_match( "forward", optarg ) )
				{
					g_opts.homotopy_param.df_dp_type = HOMOTOPY_DF_DP_FORWARD;
				}
				else if ( is_str_nocase_match( "central", optarg ) )
				{
					g_opts.homotopy_param.df_dp_type = HOMOTOPY_DF_DP_CENTRAL;
				}
				else
				{
					fprintf( stderr, "[Error] unknown ∂f/∂p approximation method %s\n", optarg );
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

