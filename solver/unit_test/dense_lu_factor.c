#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
	double x[3] = {2, 1, 3};
	int p[3];

	int n = 3;

	if ( !dense_lu_factor( n, A, p ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}
	dense_print_matrix_LU( n, A );
	dense_print_vector_i( n, p );

	dense_solve( n, A, x, p, false );
	printf( "result=\n" );
	dense_print_vector( n, x );

	return EXIT_SUCCESS;
}



