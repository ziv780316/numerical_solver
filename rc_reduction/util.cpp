//
#include <map>
#include <string>
using std::map;
using std::string;
map<string, int> g_node_id_map;

extern "C"
{
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>

#include "util.h"

static int g_n_node = 0;
};


extern "C" {

void init_node_map ()
{
	g_n_node = 0;
	g_node_id_map.clear();
	g_node_id_map["0"] = 0;
	g_node_id_map["gnd"] = 0;
}

bool is_node_exist ( char *node )
{
	string key( node );
	map<string, int>::iterator iter;
	iter = g_node_id_map.find( key );
	if ( iter == g_node_id_map.end() )
	{
		return false;
	}
	else
	{
		return true;
	}
}

int get_node_id ( char *node )
{
	string key( node );
	map<string, int>::iterator iter;
	iter = g_node_id_map.find( key );
	if ( iter == g_node_id_map.end() )
	{
		fprintf( stderr, "[ERROR] node %s doesn't exist\n", node );
		exit(1);
	}
	else
	{
		int id = g_node_id_map[key];
		return id;
	}
}

int create_node_id ( char *node )
{
	if ( is_node_exist( node ) )
	{
		fprintf( stderr, "[ERROR] node %s already exist\n", node );
		exit(1);
	}

	++g_n_node;
	int id = g_n_node;
	string key( node );
	g_node_id_map[key] = id;

	return id;
}

int get_node_map_size ()
{
	return g_n_node;
}

void dump_node_map ()
{
	for ( map<string, int>::iterator iter = g_node_id_map.begin(); iter != g_node_id_map.end(); ++iter )
	{
		printf( "%s (%d)\n", (iter->first).c_str(), iter->second );
	}
}


};
