#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	double A[6] = {1, 3, 5, 2, 4, 6}; // column-major in Fortran, row based is [1 2; 3 4; 5 6]
	double B[4] = {2, 1, 3, 4}; // column-major in Fortran, row based is [2 3; 1 4]
	double C[6] = {0.1, 0.3, 0.5, 0.2, 0.4, 0.6};

	bool a_transpose;
	bool b_transpose;
	int m = 3;
	int n = 2;
	int k = 2;
	double alpha;
	double beta;

	a_transpose = false;
	b_transpose = false;
	alpha = 1.5;
	beta = 1.0;
	dense_matrix_matrix_multiply ( m, n, k, alpha, A, B, beta, C, a_transpose, b_transpose );
	printf( "1.5*A*B\n" );
	for ( int i = 0; i < m; ++i )
	{
		for ( int j = 0; j < n; ++j )
		{
			printf( "%.10e ", *(C + j*m + i) );
		}
		printf( "\n" );
	}

	return EXIT_SUCCESS;
}

