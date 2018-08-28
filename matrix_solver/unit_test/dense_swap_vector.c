#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	double x[3] = {1, 2, 3};
	double y[3] = {0.5, 0.3, 0.1};

	int n = 3;

	printf( "before swap\n" );
	dense_print_vector( n, x );
	dense_print_vector( n, y );
	dense_swap_vector( n, x, y );
	printf( "\nafter swap\n" );
	dense_print_vector( n, x );
	dense_print_vector( n, y );

	return EXIT_SUCCESS;
}
