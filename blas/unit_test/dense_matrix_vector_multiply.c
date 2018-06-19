#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	double x[2] = {2, 1};
	double y[3] = {0.5, 0.3, 0.1};
	double A[6] = {1, 3, 5, 2, 4, 6}; // column-major in Fortran [1 2; 3 4; 5 6]

	bool transpose = false;
	int m = 3;
	int n = 2;
	double alpha = 1.0;
	double beta = 1.0;

	dense_matrix_vector_multiply ( m, n, alpha, A, x, beta, y, transpose );

	for ( int i = 0; i < m; ++i )
	{
		printf( "y%d=%.10e\n", i, y[i] );
	}

	return EXIT_SUCCESS;
}
