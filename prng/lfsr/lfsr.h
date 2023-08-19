#ifndef LFSR_H
#define LFSR_H

#include <stdint.h>
#include <stdbool.h>

typedef enum {
	LFSR_TYPE_FIBONACCI_XOR,
	LFSR_TYPE_GALOIS,
} lfsr_type_e;

typedef struct
{
	int n_taps;
	int *tap_pos; // 1-base index, record from MSB to LSB
	uint16_t *shift_expr;
} taps_t;

typedef struct
{
	taps_t taps;

	lfsr_type_e type;

	int m; // m of m-sequence
	int possible_longest_period; // 2^m - 1

	bool debug;
	bool sanity_check;

	uint16_t seed;
} lfsr_t;

void lfsr_check_illegal_state( lfsr_t *lfsr, uint16_t state );
lfsr_t *init_lfsr ( char *taps_str, char *seed_str, lfsr_type_e type, bool debug, bool sanity_check );
uint16_t lfsr_get_period( lfsr_t *lfsr );
uint16_t lfsr_get_next_state( lfsr_t *lfsr, uint16_t state );
uint16_t lfsr_str2state( lfsr_t *lfsr, char *str );
char *lfsr_state2str( lfsr_t *lfsr, uint16_t state ); // this will allocate new memory to store bits

#endif

