#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>

extern void init_node_map ();
extern bool is_node_exist ( char *node );
extern int get_node_id ( char *node );
extern int create_node_id ( char *node );
extern int get_node_map_size ();
extern void dump_node_map ();

#endif

