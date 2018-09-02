#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = REAL_NUMBER;
	double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
	double x[3] = {2, 1, 3};
	int p[3];
	int n = 3;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );
	printf( "\nx=\n" );
	dense_print_vector( n, x, type );

	if ( !dense_lu_factor( n, A, p, type ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}
	printf( "\nLU=\n" );
	dense_print_matrix_LU( n, A, type );
	printf( "\nperm=\n" );
	dense_print_vector_i( n, p, type );

	dense_solve( n, A, x, p, false, type );
	printf( "\nresult=\n" );
	dense_print_vector( n, x, type );

	return EXIT_SUCCESS;
}



