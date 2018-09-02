#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = REAL_NUMBER;
	double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9}; 
	double alpha;
	int n = 3;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );

	printf( "\nA + 0.5*I=\n" );
	alpha = 0.5;
	dense_diagonal_addition( n, A, &alpha, type );
	dense_print_matrix( n, n, A, type );

	printf( "\n2*(A + 0.5*I)=\n" );
	alpha = 2.0;
	dense_matrix_scale( n, n, A, &alpha, type );
	dense_print_matrix( n, n, A, type );

	return EXIT_SUCCESS;
}



