#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = REAL_NUMBER;
	double x[2] = {2, 1};
	double y[3] = {0.5, 0.3, 0.1};
	double y_buf[3] = {0.5, 0.3, 0.1};
	double A[6] = {1, 3, 5, 2, 4, 6}; // column-major in Fortran, row based is[1 2; 3 4; 5 6]

	bool transpose;
	int m = 3;
	int n = 2;
	double alpha;
	double beta;

	printf( "A=\n" );
	dense_print_matrix ( m, n, A, type );

	printf( "x=\n" );
	dense_print_vector ( n, x, type );

	transpose = false;
	alpha = 1.0;
	beta = 0.0;
	dense_matrix_vector_multiply ( m, n, &alpha, A, x, &beta, y, transpose, type );
	printf( "\nA*x=\n" );
	dense_print_vector ( m, y, type );

	transpose = true;
	alpha = 1.0;
	beta = 0.0;
	dense_matrix_vector_multiply ( m, n, &alpha, A, x, &beta, y, transpose, type );
	printf( "\n(A^T)*x=\n" );
	dense_print_vector ( m, y, type );

	memcpy( y, y_buf, sizeof(double) * m );
	transpose = false;
	alpha = 1.3;
	beta = 1.5;
	dense_matrix_vector_multiply ( m, n, &alpha, A, x, &beta, y, transpose, type );
	printf( "\n1.3*(A^T)*x + 1.5*y\n" );
	dense_print_vector ( m, y, type );

	return EXIT_SUCCESS;
}
