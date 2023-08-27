//
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#include "rcr.h"
#include "util.h"

static double evaluate_expr ( char *expr );
static double trigger_elapsed_timer ();
static void rcr_build_graph ( rcr_t *rcr );
static int create_node_id ( rcr_t *rcr, char *node, bool is_branch );

rcr_t *rcr_parse_netlist_spice3 ( char *netlist_in )
{
	FILE *fin = fopen( netlist_in, "r" );
	if ( !fin )
	{
		fprintf( stderr, "[ERROR] call fopen failure -> %s\n", strerror(errno) );
		exit(1);
	}
	rcr_t *rcr = (rcr_t *) calloc ( 1, sizeof(rcr_t) );
	rcr->netlist_in = strdup( netlist_in ); 
	init_node_map();

	// parse netlist
	trigger_elapsed_timer();
	printf ( "[RCR] parse netlist ...\n" );
	size_t lineno = 0;
	char buf[BUFSIZ];
	size_t len;
	char inst_name[BUFSIZ];
	char node_p_name[BUFSIZ];
	char node_n_name[BUFSIZ];
	char node_branch_name[BUFSIZ];
	char expr[BUFSIZ];
	while ( fgets(buf, BUFSIZ, fin) )
	{
		++lineno;
		len = strlen( buf );
		buf[len] = '\0';
		--len;

		// strip leading space
		char *line = NULL;
		{
			int i = 0;
			for ( ; buf[i]; ++i )
			{
				if ( !isblank(buf[i]) )
				{
					break;
				}
			}
			line = &(buf[i]);
		}

		// strip last space
		{
			int i = len - 1;
			for ( ; i >= 0; --i )
			{
				if ( !isblank(buf[i]) )
				{
					break;
				}
			}
			buf[i + 1] = '\0';
		}

		char type_ch = tolower(line[0]);
		if ( 'r' == type_ch )
		{
			int n_scanf = sscanf( line, "%s %s %s %s", inst_name, node_p_name, node_n_name, expr );
			if ( 4 != n_scanf )
			{
				fprintf( stderr, "[ERROR] r line parse failure -> line_%lu = '%s'\n", lineno, line );
				exit(1);
			}

			++(rcr->n_res);
			rcr->res_vec = (res_t *) realloc ( (void *)rcr->res_vec, rcr->n_res * sizeof(res_t) );

			res_t inst;
			inst.name = strdup( inst_name );
			inst.node_p.name = strdup( node_p_name );
			inst.node_n.name = strdup( node_n_name );
			inst.expr = strdup( expr );
			inst.r = evaluate_expr( expr );
			inst.node_p.id = create_node_id( rcr, inst.node_p.name, false );
			inst.node_n.id = create_node_id( rcr, inst.node_n.name, false );

			rcr->res_vec[rcr->n_res - 1] = inst;
		}
		else if ( 'c' == type_ch )
		{
			int n_scanf = sscanf( line, "%s %s %s %s", inst_name, node_p_name, node_n_name, expr );
			if ( 4 != n_scanf )
			{
				fprintf( stderr, "[ERROR] c line parse failure -> line_%lu = '%s'\n", lineno, line );
				exit(1);
			}

			++(rcr->n_cap);
			rcr->cap_vec = (cap_t *) realloc ( (void *)rcr->cap_vec, rcr->n_cap * sizeof(cap_t) );

			cap_t inst;
			inst.name = strdup( inst_name );
			inst.node_p.name = strdup( node_p_name );
			inst.node_n.name = strdup( node_n_name );
			inst.expr = strdup( expr );
			inst.c = evaluate_expr( expr );
			inst.node_p.id = create_node_id( rcr, inst.node_p.name, false );
			inst.node_n.id = create_node_id( rcr, inst.node_n.name, false );

			rcr->cap_vec[rcr->n_cap - 1] = inst;
		}
		else if ( 'i' == type_ch )
		{
			int n_scanf = sscanf( line, "%s %s %s %s", inst_name, node_p_name, node_n_name, expr );
			if ( 4 != n_scanf )
			{
				fprintf( stderr, "[ERROR] %c line parse failure -> line_%lu = '%s'\n", type_ch, lineno, line );
				exit(1);
			}

			++(rcr->n_src);
			++(rcr->n_isrc);
			rcr->src_vec = (src_t *) realloc ( (void *)rcr->src_vec, rcr->n_src * sizeof(src_t) );

			src_t inst;
			inst.type = SRC_TYPE_I;
			inst.name = strdup( inst_name );
			inst.node_p.name = strdup( node_p_name );
			inst.node_n.name = strdup( node_n_name );
			inst.expr = strdup( expr );
			inst.dc = evaluate_expr( expr );
			inst.node_p.id = create_node_id( rcr, inst.node_p.name, false );
			inst.node_n.id = create_node_id( rcr, inst.node_n.name, false );
			record_non_r_node( inst_name, node_p_name );
			record_non_r_node( inst_name, node_n_name );

			rcr->src_vec[rcr->n_src - 1] = inst;
		}
	}

	// parse vsrc in order to create branch in last order
	rewind( fin );
	while ( fgets(buf, BUFSIZ, fin) )
	{
		++lineno;
		len = strlen( buf );
		buf[len] = '\0';
		--len;

		// strip leading space
		char *line = NULL;
		{
			int i = 0;
			for ( ; buf[i]; ++i )
			{
				if ( !isblank(buf[i]) )
				{
					break;
				}
			}
			line = &(buf[i]);
		}

		// strip last space
		{
			int i = len - 1;
			for ( ; i >= 0; --i )
			{
				if ( !isblank(buf[i]) )
				{
					break;
				}
			}
			buf[i + 1] = '\0';
		}

		char type_ch = tolower(line[0]);
		if ( 'v' == type_ch )
		{
			int n_scanf = sscanf( line, "%s %s %s %s", inst_name, node_p_name, node_n_name, expr );
			if ( 4 != n_scanf )
			{
				fprintf( stderr, "[ERROR] %c line parse failure -> line_%lu = '%s'\n", type_ch, lineno, line );
				exit(1);
			}

			++(rcr->n_src);
			++(rcr->n_vsrc);
			rcr->src_vec = (src_t *) realloc ( (void *)rcr->src_vec, rcr->n_src * sizeof(src_t) );

			src_t inst;
			inst.type = SRC_TYPE_V;
			inst.name = strdup( inst_name );
			inst.node_p.name = strdup( node_p_name );
			inst.node_n.name = strdup( node_n_name );
			inst.expr = strdup( expr );
			inst.dc = evaluate_expr( expr );
			inst.node_p.id = create_node_id( rcr, inst.node_p.name, false );
			inst.node_n.id = create_node_id( rcr, inst.node_n.name, false );

			sprintf( node_branch_name, "%s:p", inst_name );
			inst.node_branch.name = strdup( node_branch_name );
			inst.node_branch.id = create_node_id( rcr, inst.node_branch.name, true );

			record_non_r_node( inst_name, node_p_name );
			record_non_r_node( inst_name, node_n_name );

			rcr->src_vec[rcr->n_src - 1] = inst;
		}
	}

	double elapsed_time = trigger_elapsed_timer();
	printf ( "[RCR] parse netlist spent %.3lf sec\n", elapsed_time );

	rcr_build_graph( rcr );

	return rcr;
}

void rcr_dump_db ( rcr_t *rcr )
{
	printf( "[RCR] following shows %d nodes\n", rcr->n_node );
	for ( int i = 0; i < rcr->n_node; ++i )
	{
		printf( "%s (%d) type=%c\n", rcr->node_vec[i].name, rcr->node_vec[i].id, (rcr->node_vec[i].is_branch ? 'I' : 'V')  );
	}

	printf( "[RCR] following shows %d res\n", rcr->n_res );
	for ( int i = 0; i < rcr->n_res; ++i )
	{
		res_t *inst = &(rcr->res_vec[i]);
		printf( "%d: %s %s(%d) %s(%d) %.15le\n", i + 1, inst->name, inst->node_p.name, inst->node_p.id, inst->node_n.name, inst->node_n.id, inst->r );
	}

	printf( "[RCR] following shows %d vsrc\n", rcr->n_vsrc );
	for ( int i = 0; i < rcr->n_src; ++i )
	{
		src_t *inst = &(rcr->src_vec[i]);
		if ( inst->type == SRC_TYPE_V )
		{
			printf( "%d: %s %s(%d) %s(%d) %.15le %s(%d)\n", i + 1, inst->name, inst->node_p.name, inst->node_p.id, inst->node_n.name, inst->node_n.id, inst->dc, inst->node_branch.name, inst->node_branch.id );
		}
	}

	printf( "[RCR] following shows %d r node\n", rcr->n_r_node );
	dump_r_node_map();

	printf( "[RCR] following shows %d non-r node\n", rcr->n_non_r_node );
	dump_non_r_node_map();
}

void rcr_r_reduction ( rcr_t *rcr )
{
}

static double evaluate_expr ( char *expr )
{
	return atof( expr );
}

__attribute__((constructor))
static double trigger_elapsed_timer ()
{
	static bool init = false;
	static struct timespec ts;
	if ( !init )
	{
		init = true;

		if ( -1 == clock_gettime( CLOCK_REALTIME, &ts ) )
		{
			fprintf( stderr, "[ERROR] call clock_gettime failure -> %s\n", strerror(errno) );
			exit(1);
		}

		return NAN;
	}

	static struct timespec ts_now;
	if ( -1 == clock_gettime( CLOCK_REALTIME, &ts_now ) )
	{
		fprintf( stderr, "[ERROR] call clock_gettime failure -> %s\n", strerror(errno) );
		exit(1);
	}

	double elapsed_time = (ts_now.tv_sec - ts.tv_sec) + 1e-9*(ts_now.tv_nsec - ts.tv_nsec);
	ts = ts_now;

	return elapsed_time;
}

static void rcr_build_graph ( rcr_t *rcr )
{
	trigger_elapsed_timer();
	printf ( "[RCR] build graph ...\n" );

	for ( int i = 0; i < rcr->n_res; ++i )
	{
		record_r( &(rcr->res_vec[i]) );
	}
	rcr->n_r_node = get_r_node_size();
	rcr->n_non_r_node = get_non_r_node_size();


	double elapsed_time = trigger_elapsed_timer();
	printf ( "[RCR] build graph spent %.3lf sec\n", elapsed_time );
}

static int create_node_id ( rcr_t *rcr, char *node_name, bool is_branch )
{
	if ( is_gnd( node_name ) )
	{
		return 0;
	}

	if ( is_node_exist( node_name ) )
	{
		return get_node_id( node_name );
	}

	if ( 0 == rcr->n_node )
	{
		rcr->node_vec = (node_t *) calloc ( 1, sizeof(node_t) );
		rcr->node_vec[0].name = "0";
		rcr->node_vec[0].id = 0;
		rcr->node_vec[0].is_branch = false;
		++rcr->n_node;
	}

	int id = rcr->n_node;
	insert_node_id( node_name, id );
	rcr->node_vec = (node_t *) realloc ( (void *)rcr->node_vec, (rcr->n_node + 1) * sizeof(node_t) );
	rcr->node_vec[rcr->n_node].name = node_name;
	rcr->node_vec[rcr->n_node].id = id;
	rcr->node_vec[rcr->n_node].is_branch = is_branch;
	++rcr->n_node;

	return id;
}

bool is_gnd ( char *node )
{
	if ( (0 == strcasecmp( node, "0" )) || (0 == strcasecmp(node, "gnd")) )
	{
		return true;
	}
	else
	{
		return false;
	}
}
