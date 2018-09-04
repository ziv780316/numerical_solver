#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double x[6] = {1, 2, 3, 4, 5, 6};
	double y[6] = {2, -1, 2, -5, 3, 1};
	int n = 3;
	complex_t result;

	printf( "x=\n" );
	dense_print_vector( n, x, type );

	printf( "\ny=\n" );
	dense_print_vector( n, y, type );

	bool conjugate = false;
	printf( "\nxT*y= " );
	dense_vector_inner_product( n, x, y, (double *)&result, conjugate, type );
	printf( "%.10e + i*%.10e\n", result.real, result.imag );

	conjugate = true;
	printf( "\nxH*y= " );
	dense_vector_inner_product( n, x, y, (double *)&result, conjugate, type );
	printf( "%.10e + i*%.10e\n", result.real, result.imag );

	return EXIT_SUCCESS;
}

