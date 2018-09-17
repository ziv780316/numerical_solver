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
	double y[6] = {1, -2, 1, 4, 1, 6};
	int p[3];
	int n = 3;
	int nrhs = 1;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );
	printf( "\nx=\n" );
	dense_print_vector( n, x, type );
	printf( "\ny=\n" );
	dense_print_vector( 2 * n, y, type );

	if ( !dense_lu_factor( n, A, p, type ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}
	printf( "\nLU=\n" );
	dense_print_matrix_LU( n, A, type );
	printf( "\np=\n" );
	dense_print_vector_i( n, p, type );

	dense_solve( n, nrhs, A, x, p, TRANS_NONE, type );
	printf( "\nsolve x=\n" );
	dense_print_vector( n, x, type );

	nrhs = 2;
	dense_solve( n, nrhs, A, y, p, TRANS_NONE, type );
	printf( "\nsolve y1=\n" );
	dense_print_vector( n, y, type );
	printf( "\nsolve y2=\n" );
	dense_print_vector( n, y + n, type );

	return EXIT_SUCCESS;
}



