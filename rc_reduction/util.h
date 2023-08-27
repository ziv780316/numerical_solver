#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>
#include "rcr.h"

extern bool is_gnd ( char *node );
extern void init_node_map ();
extern bool is_node_exist ( char *node );
extern int get_node_id ( char *node );
extern void insert_node_id ( char *node, int id );
extern int get_node_map_size ();
extern void dump_node_map ();
extern void record_r ( res_t *inst );
extern void record_non_r_node ( char *inst_name, char *node );
extern int get_r_node_size ();
extern int get_non_r_node_size ();
extern void dump_r_node_map ();
extern void dump_non_r_node_map ();

#endif

