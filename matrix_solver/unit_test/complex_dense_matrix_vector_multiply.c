#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double x[2*2] = {2, 1, 3, -4};
	double y[3*2] = {0};
	double y_buf[3*2] = {6, 5, 4, 3, 2, 1};
	double A[6*2] = {1, 2, 3, 4, 5, 6, -1, 2, -3, 4, -5, 6}; 

	transpose_type tran;
	int m = 3;
	int n = 2;
	complex_t alpha;
	complex_t beta;

	printf( "A=\n" );
	dense_print_matrix ( m, n, A, type );

	printf( "x=\n" );
	dense_print_vector ( n, x, type );

	printf( "y=\n" );
	dense_print_vector ( m, y_buf, type );

	tran = TRANS_NONE;
	alpha.real = 1.0;
	alpha.imag = 2.0;
	beta.real = 0.0;
	beta.imag = 0.0;
	dense_matrix_vector_multiply ( m, n, (double *)&alpha, A, x, (double *)&beta, y, tran, type );
	printf( "\n(1+2i)*A*x=\n" );
	dense_print_vector ( m, y, type );

	tran = TRANS_NORMAL;
	alpha.real = 1.0;
	alpha.imag = 2.0;
	beta.real = 0.0;
	beta.imag = 0.0;
	dense_matrix_vector_multiply ( m, m, (double *)&alpha, A, y_buf, (double *)&beta, y, tran, type );
	printf( "\n(1+2i)*(A^T)*y=\n" );
	dense_print_vector ( n, y, type );

	tran = TRANS_CONJUGATE;
	alpha.real = 1.0;
	alpha.imag = 2.0;
	beta.real = 0.0;
	beta.imag = 0.0;
	dense_matrix_vector_multiply ( m, m, (double *)&alpha, A, y_buf, (double *)&beta, y, tran, type );
	printf( "\n(1+2i)*(A^H)*y=\n" );
	dense_print_vector ( n, y, type );

	memcpy( y, y_buf, sizeof(double) * 2 * m );
	tran = TRANS_NONE;
	alpha.real = 1.0;
	alpha.imag = 2.0;
	beta.real = 2.0;
	beta.imag = -1.0;
	dense_matrix_vector_multiply ( m, n, (double *)&alpha, A, x, (double *)&beta, y, tran, type );
	printf( "\n(1+2i)*A*x + (2-i)*y\n" );
	dense_print_vector ( m, y, type );

	return EXIT_SUCCESS;
}
