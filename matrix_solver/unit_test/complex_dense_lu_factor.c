#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	number_type type = COMPLEX_NUMBER;
	double A[9*2] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 9, 8, 7, -6, -5, -4, -3, 2, 1}; 
	double x[3*2] = {2, 1, 3, 1, 9, -1};
	double y[6*2] = {1, -2, 1, 4, 1, 6, 6, 5, 1, 0, 2, 3};
	double z[6*2] = {1, -2, 1, 4, 1, 6, 6, 5, 1, 0, 2, 3};
	double w[6*2] = {1, -2, 1, 4, 1, 6, 6, 5, 1, 0, 2, 3};
	double t[3*2] = {0};
	double B[9*2] = {0};
	double C[9*2] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 9, 8, 7, -6, -5, -4, -3, 2, 1}; 
	double D[9*2] = {0};
	double E[9*2] = {0};
	int p[3] = {0};
	int n = 3;
	int nrhs = 1;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );
	printf( "\nx=\n" );
	dense_print_vector( n, x, type );
	printf( "\ny=\n" );
	dense_print_vector( 2 * n, y, type );

	printf( "\nA*(A**H)=\n" );
	complex_t alpha = { .real = 1.0, .imag = 0.0 };
	complex_t beta = { .real = 0.0, .imag = 0.0 };
	dense_matrix_matrix_multiply ( n, n, n, n, (double *)&alpha, A, C, (double *)&beta, B, TRANS_NONE, TRANS_CONJUGATE, type );
	dense_print_matrix( n, n, B, type );

	printf( "\nA*(A**T)=\n" );
	dense_matrix_matrix_multiply ( n, n, n, n, (double *)&alpha, A, C, (double *)&beta, D, TRANS_NONE, TRANS_NORMAL, type );
	dense_print_matrix( n, n, D, type );

	if ( !dense_lu_factor( n, A, p, FACTOR_LU_RIGHT_LOOKING, type ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}
	printf( "\nLU=\n" );
	dense_print_matrix_LU( n, A, type );
	printf( "\nperm=\n" );
	dense_print_vector_i( n, p );
	printf( "\nP=\n" );
	dense_print_matrix_perm( n, p );

	dense_solve( n, nrhs, A, x, p, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, type );
	printf( "\nsolve x=\n" );
	dense_print_vector( n, x, type );

	nrhs = 2;
	dense_solve( n, nrhs, A, y, p, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, type );
	printf( "\nsolve y1=\n" );
	dense_print_vector( n, y, type );
	printf( "\nsolve y2=\n" );
	dense_print_vector( n, y + 2*n, type );

	if ( !dense_lu_factor( n, B, NULL, FACTOR_LU_CHOLESKY, type ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}
	printf( "\nLU=\n" );
	dense_print_matrix_trig( n, B, TRIG_LOWER, type );
	nrhs = 2;
	dense_solve( n, nrhs, B, z, NULL, FACTOR_LU_CHOLESKY, TRANS_NONE, type );
	printf( "\nsolve y1=\n" );
	dense_print_vector( n, z, type );
	printf( "\nsolve y2=\n" );
	dense_print_vector( n, z + 2*n, type );

	memcpy( E, D, sizeof(double) * 9 * 2 );
	if ( !dense_lu_factor( n, D, p, FACTOR_LU_BUNCH_KAUFMAN, type ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}
	printf( "\nLU=\n" );
	dense_print_matrix_trig( n, D, TRIG_LOWER, type );
	nrhs = 2;
	dense_solve( n, nrhs, D, w, p, FACTOR_LU_BUNCH_KAUFMAN, TRANS_NONE, type );
	printf( "\nsolve w1=\n" );
	dense_print_vector( n, w, type );
	printf( "\nsolve w2=\n" );
	dense_print_vector( n, w + 2*n, type );
	printf( "\nsolve w1=\n" );
	printf( "\nD*s1=\n" );
	dense_matrix_vector_multiply ( n, n, (double *)&alpha, E, w, (double *)&beta, t, TRANS_NONE , type );
	dense_print_vector( n, t, type );
	printf( "\nD*s2=\n" );
	dense_matrix_vector_multiply ( n, n, (double *)&alpha, E, w + 2*n, (double *)&beta, t, TRANS_NONE , type );
	dense_print_vector( n, t, type );

	return EXIT_SUCCESS;
}



