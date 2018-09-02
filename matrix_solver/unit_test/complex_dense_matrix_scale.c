#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double A[9*2] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8, 7, 6, 5, 4, 3, 2, 1}; 
	double D[9*2];
	double alpha[2];
	int n = 3;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );

	printf( "\n(1-2i)*A=\n" );
	alpha[0] = 1.0;
	alpha[1] = -2.0;
	dense_matrix_scale( n, n, A, alpha, type );
	dense_print_matrix( n, n, A, type );

	printf( "\n(1-2i)*A + (0.5+0.1i)*I=\n" );
	alpha[0] = 0.5;
	alpha[1] = 0.1;
	dense_diagonal_addition( n, A, alpha, type );
	dense_print_matrix( n, n, A, type );

	printf( "\nD=\n" );
	dense_matrix_get_diagonal( n, A, D, type );
	dense_print_vector( n, D, type );

	return EXIT_SUCCESS;
}



