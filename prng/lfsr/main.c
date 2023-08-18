#include <stdio.h>
#include <stdlib.h>

#include "lfsr.h"

int main ( int argc, char **argv )
{
	// illegal seed if all zeroes in LFSR_TYPE_FIBONACCI_XOR
	// the initial state doesn't affect period calculation

	// 1 + x^2 + x^3
	{
		// bit1 -> bit16
		// LSB     MSB
		char *seed3_str = "111";
		char *taps3_str = "011";

		lfsr_t *lfsr = init_lfsr ( taps3_str, seed3_str, LFSR_TYPE_FIBONACCI_XOR, true, true );

		int period = lfsr_get_period( lfsr );
		printf( "period = %d (2^%d - 1 = %d)\n", period, lfsr->m, (1 << lfsr->m) - 1 );
	}

	// 1 + x^7 + x^10
	{
		// bit1 -> bit16
		// LSB     MSB
		char *seed10_str = "1100001111";
		char *taps10_str = "0000001001";

		lfsr_t *lfsr = init_lfsr ( taps10_str, seed10_str, LFSR_TYPE_FIBONACCI_XOR, true, true );

		int period = lfsr_get_period( lfsr );
		printf( "period = %d (2^%d - 1 = %d)\n", period, lfsr->m, (1 << lfsr->m) - 1 );
	}

	// 1 + x^4 + x^13 + x^15 + x^16
	{
		// bit1 -> bit16
		// LSB     MSB
		char *seed16_str = "1000000000000000";
		char *taps16_str = "0001000000001011";

		lfsr_t *lfsr = init_lfsr ( taps16_str, seed16_str, LFSR_TYPE_FIBONACCI_XOR, true, true );

		int period = lfsr_get_period( lfsr );
		printf( "period = %d (2^%d - 1 = %d)\n", period, lfsr->m, (1 << lfsr->m) - 1 );
	}

	return EXIT_SUCCESS;
}

