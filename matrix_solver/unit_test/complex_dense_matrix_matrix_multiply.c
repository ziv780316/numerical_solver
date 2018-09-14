#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double A[6*2] = {1, 2, 3, 4, 5, 6, -1, 2, -3, 4, -5, 6}; // 3 * 2
	int m = 3;
	int n = 3;
	int k = 2;
	double B[6*2] = {2, 1, 3, 4, 0, 6, 7, 1, 2, 3, -3, -1};  // 2 * 3
	double C[9*2] = {0}; // 3 * 3
	double D[4*2] = {1, 2, 3, 4, 5, 6, -7, -8}; // 2 * 2

	transpose_type a_transpose;
	transpose_type b_transpose;
	complex_t alpha;
	complex_t beta;

	printf( "A\n" );
	dense_print_matrix( m, k, A, type );

	printf( "\nB\n" );
	dense_print_matrix( k, n, B, type );

	printf( "\nC\n" );
	dense_print_matrix( m, n, C, type );

	printf( "\nD\n" );
	dense_print_matrix( k, k, D, type );

	a_transpose = TRANS_NONE;
	b_transpose = TRANS_NONE;
	alpha.real = 1.0;
	alpha.imag = 2.0;
	beta.real = 0.0;
	beta.imag = 0.0;
	dense_matrix_matrix_multiply ( m, k, k, n, (double *)&alpha, A, B, (double *)&beta, C, a_transpose, b_transpose, type );
	printf( "\nC = (1+2*i)*A*B\n" );
	dense_print_matrix( m, n, C, type );

	a_transpose = TRANS_CONJUGATE;
	b_transpose = TRANS_NORMAL;
	alpha.real = 1.0;
	alpha.imag = 2.0;
	beta.real = 3.0;
	beta.imag = -1.0;
	dense_matrix_matrix_multiply ( m, k, k, n, (double *)&alpha, A, B, (double *)&beta, D, a_transpose, b_transpose, type );
	printf( "\nD = (1+2*i)*(A**H)*(B**T) + (3-i)*D\n" );
	dense_print_matrix( k, k, D, type );

	return EXIT_SUCCESS;
}

