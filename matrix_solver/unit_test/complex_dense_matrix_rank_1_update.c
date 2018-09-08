#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double x[2*2] = {1, 2, 3, 4};
	double y[3*2] = {2, -1, 2, -5, 3, 1};
	double A[2*3*2] = {0};
	int m = 2;
	int n = 3;

	printf( "x=\n" );
	dense_print_vector( m, x, type );

	printf( "\ny=\n" );
	dense_print_vector( n, y, type );

	bool conjugate = false;
	printf( "\nx*yT=\n" );
	dense_vector_outer_product( m, n, x, y, A, conjugate, type );
	dense_print_matrix( m, n, A, type );

	conjugate = true;
	complex_t alpha = { .real = 1.0, .imag= 2.0 };
	printf( "\nA + x*yH=\n" );
	dense_maxtrix_rank_1_update( m, n, A, (double *)&alpha, x, y, conjugate, type );
	dense_print_matrix( m, n, A, type );

	return EXIT_SUCCESS;
}

