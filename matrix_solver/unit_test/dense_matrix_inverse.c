#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
	double B[9] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
	double C[9] = {0};
	int p[3];
	int n = 3;

	printf( "A=\n" );
	dense_print_matrix( n, n, A );

	if ( !dense_matrix_inverse( n, A, p ) )
	{
		fprintf( stderr, "[Error] inverse matrix fail\n" );
		abort();
	}

	printf( "\nA^-1=\n" );
	dense_print_matrix( n, n, A );

	bool a_transpose = false;
	bool b_transpose = false;
	double alpha = 1.0;
	double beta = 1.0;
	dense_matrix_matrix_multiply ( n, n, n, alpha, B, A, beta, C, a_transpose, b_transpose );
	printf( "\nA*A^-1=\n" );
	dense_print_matrix( n, n, C );

	return EXIT_SUCCESS;
}



