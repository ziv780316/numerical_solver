#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = REAL_NUMBER;
	double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
	double B[9] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
	double C[9] = {0};
	double D[9] = {0};
	double E[9] = {0};
	int p[3] = {0};
	int n = 3;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );

	if ( !dense_matrix_inverse( n, A, p, FACTOR_LU_RIGHT_LOOKING, type ) )
	{
		fprintf( stderr, "[Error] inverse matrix fail\n" );
		abort();
	}

	printf( "\np=\n" );
	dense_print_vector_i( n, p );

	printf( "\nA^-1=\n" );
	dense_print_matrix( n, n, A, type );

	double alpha = 1.0;
	double beta = 0.0;
	dense_matrix_matrix_multiply ( n, n, n, n, &alpha, B, A, &beta, C, TRANS_NONE, TRANS_NONE, type );
	printf( "\nA*A^-1=\n" );
	dense_print_matrix( n, n, C, type );
	
	printf( "\nA*(A**T)=\n" );
	dense_matrix_matrix_multiply ( n, n, n, n, &alpha, B, B, &beta, D, TRANS_NONE, TRANS_NORMAL, type );
	dense_print_matrix( n, n, D, type );
	if ( !dense_matrix_inverse( n, D, p, FACTOR_LU_CHOLESKY, type ) )
	{
		fprintf( stderr, "[Error] inverse matrix fail\n" );
		abort();
	}
	printf( "\n(A*A^T)^-1=\n" );
	dense_print_matrix( n, n, D, type );

	dense_matrix_matrix_multiply ( n, n, n, n, &alpha, B, B, &beta, D, TRANS_NONE, TRANS_NORMAL, type );
	if ( !dense_matrix_inverse( n, D, p, FACTOR_LU_BUNCH_KAUFMAN, type ) )
	{
		fprintf( stderr, "[Error] inverse matrix fail\n" );
		abort();
	}
	printf( "\n(A*A^T)^-1=\n" );
	dense_print_matrix( n, n, D, type );

	memcpy( E, D, sizeof(double) * 9 );
	dense_matrix_matrix_multiply ( n, n, n, n, &alpha, B, B, &beta, D, TRANS_NONE, TRANS_NORMAL, type );
	printf( "\n(A*A^T)*(A*A^T)^-1=\n" );
	dense_matrix_matrix_multiply ( n, n, n, n, &alpha, D, E, &beta, C, TRANS_NONE, TRANS_NONE, type );
	dense_print_matrix( n, n, C, type );

	return EXIT_SUCCESS;
}



