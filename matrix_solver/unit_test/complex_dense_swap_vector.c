#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double x[3*2] = {1, 2, 3, 3, 2, 1};
	double y[3*2] = {0.5, 0.3, 0.1, 1, 3, 1};

	int n = 3;

	printf( "before swap\n" );
	printf( "x=\n" );
	dense_print_vector( n, x, type );
	printf( "\ny=\n" );
	dense_print_vector( n, y, type );

	dense_swap_vector( n, x, y, type );

	printf( "\nafter swap\n" );
	printf( "x=\n" );
	dense_print_vector( n, x, type );
	printf( "\ny=\n" );
	dense_print_vector( n, y, type );

	return EXIT_SUCCESS;
}
