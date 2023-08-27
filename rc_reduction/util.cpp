//
extern "C"
{
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>

#include "util.h"
};

#include <map>
#include <string>
#include <vector>
using std::map;
using std::string;
using std::vector;
map<string, int> g_node_id_map;
map<string, vector<res_t *> > g_r_node_map;
map<string, vector<string> > g_non_r_node_map;
extern "C" {

void init_node_map ()
{
	g_node_id_map.clear();
	g_node_id_map["0"] = 0;
	g_node_id_map["gnd"] = 0;

	g_r_node_map.clear();
	g_non_r_node_map.clear();
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

int get_node_map_size ()
{
	return g_node_id_map.size();
}

void insert_node_id ( char *node, int id )
{
	string key( node );
	g_node_id_map[key] = id;
}

void dump_node_map ()
{
	for ( map<string, int>::iterator iter = g_node_id_map.begin(); iter != g_node_id_map.end(); ++iter )
	{
		printf( "%s (%d)\n", (iter->first).c_str(), iter->second );
	}
}

int get_r_node_size ()
{
	return g_r_node_map.size();
}

void record_r ( res_t *inst )
{
	// pos
	if ( !is_gnd(inst->node_p.name) )
	{
		string key( inst->node_p.name );
		g_r_node_map[key].push_back( inst );
	}

	// neg
	if ( !is_gnd(inst->node_n.name) )
	{
		string key( inst->node_n.name );
		g_r_node_map[key].push_back( inst );
	}
}

void dump_r_node_map ()
{
	int cnt = 0;
	for ( map<string, vector<res_t *> >::iterator iter = g_r_node_map.begin(); iter != g_r_node_map.end(); ++iter )
	{
		++cnt;
		printf( "%d: %s connect to %lu r\n", cnt, iter->first.c_str(), iter->second.size() );
		for ( int i = 0; i < iter->second.size(); ++i )
		{
			printf( "  %d: %s %.15le\n", i + 1, iter->second[i]->name, iter->second[i]->r );
		}
	}
}

void record_non_r_node ( char *inst_name, char *node )
{
	if ( !is_gnd(node) )
	{
		string key( node );
		string val( inst_name );
		g_non_r_node_map[key].push_back( val );
	}
}

int get_non_r_node_size ()
{
	return g_non_r_node_map.size();
}

void dump_non_r_node_map ()
{
	int cnt = 0;
	for ( map<string, vector<string> >::iterator iter = g_non_r_node_map.begin(); iter != g_non_r_node_map.end(); ++iter )
	{
		++cnt;
		printf( "%d: %s connect to %lu non-r inst\n", cnt, iter->first.c_str(), iter->second.size() );
		for ( int i = 0; i < iter->second.size(); ++i )
		{
			printf( "  %d: %s\n", i + 1, iter->second[i].c_str() );
		}
	}
}

}
