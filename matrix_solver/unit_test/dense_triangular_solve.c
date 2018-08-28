#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	double A[9] = {1, 2, 3, 0, 4, 5, 0, 0, 6}; // column-major in Fortran, row based is [1; 2 4; 3 5 6]
	double B[9] = {1, 0, 0, 2, 4, 0, 3, 5, 6}; 
	double b[3] = {3, 2, 1}; 
	double b_buf[3] = {3, 2, 1}; 

	bool is_lower_triangular;
	bool transpose;
	bool is_unit_triangular;
	int n = 3;

	is_lower_triangular = true;
	transpose = false;
	is_unit_triangular = false;
	dense_triangular_solve( n, A, b, is_lower_triangular, transpose, is_unit_triangular );
	printf( "A=\n" );
	dense_print_matrix( n, n, A );
	for ( int i = 0; i < n; ++i )
	{
		printf( "x%d = %.10e\n", i, b[i] );
	}

	memcpy( b, b_buf, sizeof(double) * n );
	is_lower_triangular = false;
	dense_triangular_solve( n, B, b, is_lower_triangular, transpose, is_unit_triangular );
	for ( int i = 0; i < n; ++i )
	{
		printf( "x%d = %.10e\n", i, b[i] );
	}

	memcpy( b, b_buf, sizeof(double) * n );
	is_lower_triangular = true;
	is_unit_triangular = true;
	dense_triangular_solve( n, A, b, is_lower_triangular, transpose, is_unit_triangular );
	for ( int i = 0; i < n; ++i )
	{
		printf( "x%d = %.10e\n", i, b[i] );
	}

	return EXIT_SUCCESS;
}


