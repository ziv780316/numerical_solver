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
	char *node_p;
	char *node_n;
	char *expr;

	int node_p_id;
	int node_n_id;

	double r;
} res_t;

typedef struct 
{
	char *name;
	char *node_p;
	char *node_n;
	char *expr;

	int node_p_id;
	int node_n_id;

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
	cap_t *cap_vec;
	res_t *res_vec;
	int n_node; // exclude gnd
	int n_res;
	int n_cap;
	int n_cpl_cap;
	int n_gnd_cap;
	int n_decpl_cap;

	// IO
	FILE *fout;
	char *netlist_in;
	char *netlist_out;
} rcr_t;

rcr_t *rcr_parse_netlist_spice3 ( char *netlist_in ); 
void rcr_dump_db ( rcr_t *rcr );


#endif

