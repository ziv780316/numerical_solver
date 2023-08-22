//
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <ctype.h>

#include "rcr.h"
#include "util.h"

static double evaluate_expr ( char *expr );
extern int create_node_id ( char *node );

rcr_t *rcr_parse_netlist_spice3 ( char *netlist_in )
{
	FILE *fin = fopen( netlist_in, "r" );
	if ( !fin )
	{
		fprintf( stderr, "[ERROR] fopen failure -> %s\n", strerror(errno) );
		exit(1);
	}
	rcr_t *rcr = (rcr_t *) calloc ( 1, sizeof(rcr_t) );
	rcr->netlist_in = strdup( netlist_in ); 
	init_node_map();

	size_t lineno = 0;
	char buf[BUFSIZ];
	size_t len;
	char name[BUFSIZ];
	char node_p[BUFSIZ];
	char node_n[BUFSIZ];
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

		if ( 'r' == tolower(line[0]) )
		{
			int n_scanf = sscanf( line, "%s %s %s %s", name, node_p, node_n, expr );
			if ( 4 != n_scanf )
			{
				fprintf( stderr, "[ERROR] r line parse failure -> line_%lu = '%s'\n", lineno, line );
				exit(1);
			}

			++(rcr->n_res);
			rcr->res_vec = (res_t *) realloc ( (void *)rcr->res_vec, rcr->n_res * sizeof(res_t) );

			res_t inst;
			inst.name = strdup( name );
			inst.node_p = strdup( node_p );
			inst.node_n = strdup( node_n );
			inst.expr = strdup( expr );
			inst.r = evaluate_expr( expr );
			if ( !is_node_exist( node_p ) )
			{
				inst.node_p_id = create_node_id( node_p );
			}
			else
			{
				inst.node_p_id = get_node_id( node_p );
			}
			if ( !is_node_exist( node_n ) )
			{
				inst.node_n_id = create_node_id( node_n );
			}
			else
			{
				inst.node_n_id = get_node_id( node_n );
			}

			rcr->res_vec[rcr->n_res - 1] = inst;
		}
	}

	rcr->n_node = get_node_map_size();
	
	return rcr;
}

void rcr_dump_db ( rcr_t *rcr )
{
	printf( "[RCR] following shows %d nodes\n", rcr->n_node );
	dump_node_map();

	printf( "[RCR] following shows %d res\n", rcr->n_res );
	for ( int i = 0; i < rcr->n_res; ++i )
	{
		res_t *inst = &(rcr->res_vec[i]);
		printf( "%d: %s %s(%d) %s(%d) %.15le\n", i + 1, inst->name, inst->node_p, inst->node_p_id, inst->node_n, inst->node_n_id, inst->r );
	}
}

static double evaluate_expr ( char *expr )
{
	return atof( expr );
}
