#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = REAL_NUMBER;
	double A[6] = {1, 3, 2, 4, 5, 7}; 
	double x[2] = {1, 2};
	double y[3] = {3, 4, 5};
	double B[6] = {0};

	printf( "x=\n" );
	dense_print_vector( 2, x, type );

	printf( "\ny=\n" );
	dense_print_vector( 3, y, type );

	printf( "x * yT=\n" );
	dense_vector_outer_product( 2, 3, x, y, B, false, type );
	dense_print_matrix( 2, 3, B, type );

	printf( "\nA=\n" );
	dense_print_matrix( 2, 3, A, type );

	printf( "\nA_rank1_update=\n" );
	double alpha = 1.0;
	dense_maxtrix_rank_1_update ( 2, 3, A, &alpha, x, y, false, type );
	dense_print_matrix( 2, 3, A, type );


	return EXIT_SUCCESS;
}

