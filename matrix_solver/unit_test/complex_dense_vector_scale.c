#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double x[6] = {1, 2, 3, 4, 5, 6};
	double x_buf[6] = {1, 2, 3, 4, 5, 6};
	double alpha[2] = {0.0, 0.0};
	int n = 3;

	printf( "x=\n" );
	dense_print_vector( n, x, type );

	printf( "\nx*2=\n" );
	alpha[0] = 2.0;
	alpha[1] = 0.0;
	dense_vector_scale( n, x, alpha, type );
	dense_print_vector( n, x, type );

	memcpy( x, x_buf, sizeof(double) * 2 * n );

	printf( "\nx*(1-2i)=\n" );
	alpha[0] = 1.0;
	alpha[1] = -2.0;
	dense_vector_scale( n, x, alpha, type );
	dense_print_vector( n, x, type );

	return EXIT_SUCCESS;
}

