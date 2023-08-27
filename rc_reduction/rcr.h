//
#ifndef RCR_H
#define RCR_H

#include <stdio.h>
#include <stdbool.h>

typedef enum
{
	RCR_REDUCE_R,
	RCR_REDUCE_RC,
	RCR_REDUCE_CPL_C,
} rcr_e;

typedef struct 
{
	char *name;
	int id;
	bool is_branch;
} node_t;

typedef struct 
{
	char *name;
	node_t node_p;
	node_t node_n;
	node_t node_branch;
	char *expr;

	int type; 
#define SRC_TYPE_V 0
#define SRC_TYPE_I 1

	double dc;
} src_t;

typedef struct 
{
	char *name;
	node_t node_p;
	node_t node_n;
	char *expr;

	double r;
} res_t;

typedef struct 
{
	char *name;
	node_t node_p;
	node_t node_n;
	char *expr;

	double c;
} cap_t;

typedef struct
{
	bool debug; 

	rcr_e reduction_type;

	// reduction options
	double cpl_min;
	double res_min;

	// RC db
	node_t *node_vec;
	cap_t *cap_vec;
	res_t *res_vec;
	src_t *src_vec;
	int n_node; // exclude gnd
	int n_res;
	int n_cap;
	int n_src;
	int n_vsrc;
	int n_isrc;
	int n_cpl_cap;
	int n_gnd_cap;
	int n_decpl_cap;
	int n_r_node;
	int n_non_r_node;
	int n_only_r_node;

	// IO
	FILE *fout;
	char *netlist_in;
	char *netlist_out;
} rcr_t;

rcr_t *rcr_parse_netlist_spice3 ( char *netlist_in ); 
void rcr_r_reduction ( rcr_t *rcr );
void rcr_dump_db ( rcr_t *rcr );
bool is_gnd ( char *node );


#endif

