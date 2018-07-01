#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9}; 
	int n = 3;

	dense_diagonal_addition( n, A, 0.5 );
	dense_print_matrix( n, n, A );
	dense_matrix_scale( n, n, A, 2.0 );
	dense_print_matrix( n, n, A );

	return EXIT_SUCCESS;
}



