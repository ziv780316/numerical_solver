#include <stdio.h>
#include <stdlib.h>

#include "lfsr.h"

int main ( int argc, char **argv )
{
	// illegal seed if all zeroes in LFSR_TYPE_FIBONACCI_XOR
	// the initial state doesn't affect period calculation

	// x^3 + x^2 + 1
	{
		// bit1 -> bit16
		// LSB     MSB
		char *seed3_str = "111";
		char *taps3_str = "110";

		lfsr_t *lfsr = init_lfsr ( taps3_str, seed3_str, LFSR_TYPE_FIBONACCI_XOR, true, true );

		int period = lfsr_get_period( lfsr );
		printf( "period = %d (2^%d - 1 = %d)\n", period, lfsr->m, (1 << lfsr->m) - 1 );
	}

	// x^4 + x^3 + 1
	{
		// bit1 -> bit16
		// LSB     MSB
		char *seed3_str = "1111";
		char *taps3_str = "1100";

		lfsr_t *lfsr = init_lfsr ( taps3_str, seed3_str, LFSR_TYPE_FIBONACCI_XOR, true, true );

		int period = lfsr_get_period( lfsr );
		printf( "period = %d (2^%d - 1 = %d)\n", period, lfsr->m, (1 << lfsr->m) - 1 );
	}

	// x^10 + x^7 + 1
	{
		// bit1 -> bit16
		// LSB     MSB
		char *seed10_str = "1111111111";
		char *taps10_str = "1001000000";

		lfsr_t *lfsr = init_lfsr ( taps10_str, seed10_str, LFSR_TYPE_FIBONACCI_XOR, true, true );

		int period = lfsr_get_period( lfsr );
		printf( "period = %d (2^%d - 1 = %d)\n", period, lfsr->m, (1 << lfsr->m) - 1 );
	}

	// x^16 + x^15 + x^13 + x^4 + 1
	{
		// bit1 -> bit16
		// LSB     MSB
		char *seed16_str = "1111111111111111";
		char *taps16_str = "1101000000001000";

		lfsr_t *lfsr = init_lfsr ( taps16_str, seed16_str, LFSR_TYPE_FIBONACCI_XOR, true, true );

		int period = lfsr_get_period( lfsr );
		printf( "period = %d (2^%d - 1 = %d)\n", period, lfsr->m, (1 << lfsr->m) - 1 );
	}

	return EXIT_SUCCESS;
}

