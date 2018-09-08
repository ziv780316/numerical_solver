#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double x[6] = {1, 2, 3, 4, 5, 6};
	double y[3] = {1, 2, 3};
	int n = 3;
	double norm;

	printf( "x=\n" );
	dense_print_vector( n, x, type );

	printf( "\nx_norm_max = " );
	dense_vector_norm( -1, n, x, &norm, type );
	printf( "%.10e\n", norm );
	for ( int i = 1; i <= 3; ++i )
	{
		printf( "\nx_norm_%d = ", i );
		dense_vector_norm( i, n, x, &norm, type );
		printf( "%.10e\n", norm );
	}

	type = REAL_NUMBER;
	printf( "\ny=\n" );
	dense_print_vector( n, y, type );

	printf( "\ny_norm_max = " );
	dense_vector_norm( -1, n, y, &norm, type );
	printf( "%.10e\n", norm );
	for ( int i = 1; i <= 3; ++i )
	{
		printf( "\ny_norm_%d = ", i );
		dense_vector_norm( i, n, y, &norm, type );
		printf( "%.10e\n", norm );
	}


	return EXIT_SUCCESS;
}

