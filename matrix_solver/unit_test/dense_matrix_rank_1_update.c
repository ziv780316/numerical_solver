#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = REAL_NUMBER;
	double A[4] = {1, 3, 2, 4}; 
	double x[2] = {1, 2};
	double y[2] = {3, 4};
	double B[4] = {0};

	printf( "x=\n" );
	dense_print_vector( 2, x, type );

	printf( "\ny=\n" );
	dense_print_vector( 2, y, type );

	printf( "x * yT=\n" );
	dense_vector_outer_product( 2, x, y, B, type );
	dense_print_matrix( 2, 2, B, type );

	printf( "\nA=\n" );
	dense_print_matrix( 2, 2, A, type );

	printf( "\nA_rank1_update=\n" );
	double alpha = 1.0;
	dense_maxtrix_rank_1_update ( 2, A, &alpha, x, y, type );
	dense_print_matrix( 2, 2, A, type );


	return EXIT_SUCCESS;
}

