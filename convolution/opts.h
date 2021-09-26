#ifndef OPTS_H
#define OPTS_H

#include <stdbool.h>

typedef struct
{
	bool debug;

	char *freq_response_file;
	char *input_signal_time_domain_file;
	char *output_signal_time_domain_file;
	char *configure_file;

} opt_t;

extern void show_help ();
extern void parse_cmd_options ( int argc, char **argv );
extern opt_t g_opts;

#endif
