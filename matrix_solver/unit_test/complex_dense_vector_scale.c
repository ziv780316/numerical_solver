#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	double x[6] = {1, 2, 3, 4, 5, 6};
	double x_buf[6] = {1, 2, 3, 4, 5, 6};
	double alpha = 2.0;
	int n = 3;

	printf( "x=\n" );
	complex_dense_print_vector( n, x );

	complex_dense_vector_scale( n, x, alpha );

	printf( "\nx*2=\n" );
	complex_dense_print_vector( n, x );

	memcpy( x, x_buf, sizeof(double) * 2 * n );

	printf( "\nx*(1-2i)=\n" );
	complex_dense_vector_scale_complex( n, x, 1, -2 );
	complex_dense_print_vector( n, x );

	return EXIT_SUCCESS;
}

