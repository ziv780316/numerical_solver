#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lfsr.h"

// convert state to string format
// str[15] -> str[0]
// bit16 -> bit1
// MSB     LSB
char *lfsr_state2str ( lfsr_t *lfsr, uint16_t state )
{
	char *str = (char *) calloc( lfsr->m + 1, 1 ); // + 1 is to store'\0'

	int bit;
	for ( int i = 0; i < lfsr->m; ++i )
	{
		bit = (state >> i) & 0x1;
		str[lfsr->m - 1 - i] = bit + '0';
	}

	if ( lfsr->sanity_check )
	{
		lfsr->sanity_check = false;

		uint16_t state_check = lfsr_str2state( lfsr, str );
		if ( state_check != state )
		{
			fprintf( stderr, "[ERROR] state=%0#6hx != lfsr_str2state(%s)=%0#6hx\n", state, str, state_check );
			exit(1);
		}

		lfsr->sanity_check = true;
	}

	return str;
}

// convert state to string format
// str[15] -> str[0]
// bit16 -> bit1
// MSB     LSB
uint16_t lfsr_str2state ( lfsr_t *lfsr, char *str )
{
	uint16_t state = 0; 

	// convert bit string to uint16_t
	for ( int i = 0; str[i]; ++i )
	{
		if ( '1' == str[i] )
		{
			state |= (1 << lfsr->m - 1 - i);
		}
	}

	if ( lfsr->sanity_check )
	{
		lfsr->sanity_check = false;

		char *str_check = lfsr_state2str( lfsr, state );
		if ( 0 != strncmp( str, str_check, lfsr->m ) )
		{
			fprintf( stderr, "[ERROR] str=%s != lfsr_state2str(%0#6hx)=%s\n", str, state, str_check );
			exit(1);
		}
		free( str_check );

		lfsr->sanity_check = true;
	}

	return state;
}

lfsr_t *init_lfsr ( char *taps_str, char *seed_str, lfsr_type_e type, bool debug, bool sanity_check )
{
	lfsr_t *lfsr = (lfsr_t *) calloc ( 1, sizeof(lfsr_t) );
	lfsr->type = type;
	lfsr->debug = debug;
	lfsr->sanity_check = sanity_check;

	// assert first tap bit must be MSB 
	if ( '1' != taps_str[0] )
	{
		fprintf( stderr, "[ERROR] taps_str[0] should be '1'\n" );
		exit(1);
	}
	lfsr->m = strlen( taps_str );
	lfsr->possible_longest_period = (1 << lfsr->m) - 1;

	// parse taps
	taps_t taps = { .n_taps = 0 };
	for ( int i = 0; taps_str[i]; ++i )
	{
		if ( '1' == taps_str[i] )
		{
			++(taps.n_taps);
		}
	}

	if ( 0 == taps.n_taps )
	{
		fprintf( stderr, "[ERROR] n_taps=%d should > 0\n", taps.n_taps );
		exit(1);
	}

	taps.tap_pos = (int *) calloc ( taps.n_taps, sizeof(int) );
	taps.shift_expr = (uint16_t *) calloc ( taps.n_taps, sizeof(uint16_t) );

	// record tap from MSB to LSB
	int cnt = 0;
	for ( int i = 0; taps_str[i]; ++i )
	{
		if ( '1' == taps_str[i] )
		{
			taps.tap_pos[cnt] = lfsr->m - i;;
			++cnt;
		}
	}

	lfsr->taps = taps;

	// parse seed
	lfsr->seed = lfsr_str2state( lfsr, seed_str );

	// check start state is legal or not
	lfsr_check_illegal_state( lfsr, lfsr->seed );

	// check initialize results
	if ( lfsr->debug )
	{
		printf( "---------------------------------------------\n" );
		printf( "following shows LFSR initialization results:\n" );
		printf( "seed = %.*s (%0#6hx)\n", lfsr->m, seed_str, lfsr->seed );
		printf( "taps = %.*s (m=%d)\n", lfsr->m, taps_str, lfsr->m );
		printf( "possible_longest_period = %d\n", lfsr->possible_longest_period );
		printf( "---------------------------------------------\n" );
	}

	return lfsr;
}

void lfsr_check_illegal_state( lfsr_t *lfsr, uint16_t state )
{
	if ( LFSR_TYPE_FIBONACCI_XOR == lfsr->type )
	{
		if ( 0 == state )
		{
			fprintf( stderr, "[ERROR] LFSR state cannot be all zeros\n" );
			exit(1);
		}
	}
}

uint16_t lfsr_get_period( lfsr_t *lfsr )
{
	// feedback shift process
	int period = 0;
	uint16_t state = lfsr->seed;
	while ( 1 )
	{
		++period;

		uint16_t state_next = lfsr_get_next_state( lfsr, state );
		if ( state_next == lfsr->seed )
		{
			break; 
		}
		state = state_next;

		if ( period > lfsr->possible_longest_period )
		{
			fprintf( stderr, "[ERROR] couldn't repeat state within the possible longest period %d (2^%d -1)\n", lfsr->possible_longest_period, lfsr->m );
			exit(1);
		}
	};

	return period;
}

uint16_t lfsr_get_next_state( lfsr_t *lfsr, uint16_t state )
{
	uint16_t state_next = 0;

	if ( LFSR_TYPE_FIBONACCI_XOR == lfsr->type )
	{
		int output_bit = 0;
		int shift;

		// evaluate shift operatorfrom MSB to LSB
		for ( int i = 0; i < (lfsr->taps).n_taps; ++i )
		{
			shift = (lfsr->m - (lfsr->taps).tap_pos[i]);
			(lfsr->taps).shift_expr[i] = state >> shift ;
		}

		int xor_expr = (lfsr->taps).shift_expr[0];
		for ( int i = 1; i < (lfsr->taps).n_taps; ++i )
		{
			xor_expr ^= (lfsr->taps).shift_expr[i];
		}
		output_bit = xor_expr & 0x1; // mod 2

		state_next = (state >> 1) | (output_bit << (lfsr->m - 1));

		if ( lfsr->debug )
		{
			char *bit_str;
			bit_str = lfsr_state2str( lfsr, state );
			printf( "state = %s (%0#6hx)\n", bit_str, state );
			free( bit_str );

			for ( int i = 0; i < (lfsr->taps).n_taps; ++i )
			{
				shift = (lfsr->m - (lfsr->taps).tap_pos[i]);
				bit_str = lfsr_state2str( lfsr, (lfsr->taps).shift_expr[i] );
				printf( "state >> %d = %s (%0#6hx)\n", shift, bit_str, (lfsr->taps).shift_expr[i] );
				free( bit_str );
			}

			bit_str = lfsr_state2str( lfsr, xor_expr );
			printf( "xor_expr = %s (%0#6hx)\n", bit_str, xor_expr );
			free( bit_str );

			printf( "output_bit = %d\n", output_bit );

			bit_str = lfsr_state2str( lfsr, state_next );
			printf( "state_next = %s (%0#6hx)\n", bit_str, state_next );
			free( bit_str );
		}
	}
	else if ( LFSR_TYPE_GALOIS == lfsr->type )
	{
		fprintf( stderr, "[ERROR] lfsr_type=Galois hasn't implement yet\n" );
		exit(1);
	}
	else
	{
		fprintf( stderr, "[ERROR] unknown lfsr_type=%d\n", lfsr->type );
		exit(1);
	}

	return state_next;
}

