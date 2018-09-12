#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = REAL_NUMBER;
	double A[6] = {1, 3, 5, 2, 4, 6}; // column-major in Fortran, row based is [1 2; 3 4; 5 6]
	double B[4] = {2, 1, 3, 4}; // column-major in Fortran, row based is [2 3; 1 4]
	double C[6] = {0.1, 0.3, 0.5, 0.2, 0.4, 0.6};
	double D[6] = {1, 2, 3, 4, 5, 6};
	double E[4] = {1, 2, 3, 4};
	double I[4] = {1, 0, 0, 1};

	transpose_type a_transpose;
	transpose_type b_transpose;
	double alpha;
	double beta;

	printf( "A=\n" );
	dense_print_matrix( 3, 2, A, type );

	printf( "\nB=\n" );
	dense_print_matrix( 2, 2, B, type );

	printf( "\nC=\n" );
	dense_print_matrix( 3, 2, C, type );

	printf( "\nD=\n" );
	dense_print_matrix( 3, 2, D, type );

	printf( "\nD2=\n" );
	dense_print_matrix( 2, 3, D, type );

	printf( "\nE=\n" );
	dense_print_matrix( 2, 2, E, type );

	a_transpose = TRANS_NONE;
	b_transpose = TRANS_NONE;
	alpha = 1.5;
	beta = 1.0;
	dense_matrix_matrix_multiply ( 3, 2, 2, 2, &alpha, A, B, &beta, C, a_transpose, b_transpose, type );
	printf( "\nC = 1.5*A*B + C\n" );
	dense_print_matrix( 3, 2, C, type );

	a_transpose = TRANS_NORMAL;
	b_transpose = TRANS_NONE;
	alpha = 1.0;
	beta = 1.0;
	dense_matrix_matrix_multiply ( 3, 2, 3, 2, &alpha, A, D, &beta, E, a_transpose, b_transpose, type );
	printf( "\nE = (A**T)*D + E\n" );
	dense_print_matrix( 2, 2, E, type );

	a_transpose = TRANS_NORMAL;
	b_transpose = TRANS_NORMAL;
	alpha = 1.0;
	beta = 1.0;
	dense_matrix_matrix_multiply ( 3, 2, 2, 3, &alpha, A, D, &beta, I, a_transpose, b_transpose, type );
	printf( "\nE = (A**T)*(D2**T) + I\n" );
	dense_print_matrix( 2, 2, I, type );

	return EXIT_SUCCESS;
}

