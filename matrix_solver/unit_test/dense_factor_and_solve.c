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
	int p[3] = {0};
	int n = 3;
	int nrhs = 1;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );
	printf( "\nx=\n" );
	dense_print_vector( n, x, type );

	if ( !dense_factor_and_solve( n, nrhs, A, x, p, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, type ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}

	printf( "\nLU=\n" );
	dense_print_matrix_LU( n, A, type );
	printf( "\nsolve result=\n" );
	dense_print_vector( n, x, type );

	return EXIT_SUCCESS;
}



