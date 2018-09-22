#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double A[9*2] = {1, 2, 3, 4, 5, 6, 7, 8, 10, -1, -2, -5, -6, -6, -9, 3, 2, 1}; 
	double B[9*2] = {1, 2, 3, 4, 5, 6, 7, 8, 10, -1, -2, -5, -6, -6, -9, 3, 2, 1}; 
	double C[9*2] = {0};
	int p[3] = {0};
	int n = 3;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );

	printf( "\np=\n" );
	dense_print_vector_i( n, p, type );

	if ( !dense_matrix_inverse( n, A, p, FACTOR_LU_RIGHT_LOOKING, type ) )
	{
		fprintf( stderr, "[Error] inverse matrix fail\n" );
		abort();
	}

	printf( "\nA^-1=\n" );
	dense_print_matrix( n, n, A, type );

	double alpha = 1.0;
	double beta = 0.0;
	dense_matrix_matrix_multiply ( n, n, n, n, &alpha, B, A, &beta, C, TRANS_NONE, TRANS_NONE, type );
	printf( "\nA*A^-1=\n" );
	dense_print_matrix( n, n, C, type );

	return EXIT_SUCCESS;
}



