//
#include <stdio.h>
#include <stdlib.h>

#include "rcr.h"

int main ( int argc, char **argv )
{
	rcr_t *rcr = rcr_parse_netlist_spice3( "input.sp" );
	rcr_dump_db ( rcr );

	return EXIT_SUCCESS;
}

