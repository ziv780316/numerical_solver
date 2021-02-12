#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	g_matrix_print_format = MATRIX_PRINT_FORMAT_MATLAB;
	number_type type = REAL_NUMBER;
	double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 10}; 
	double x[3] = {2, 1, 3};
	double y[6] = {1, -2, 1, 4, 1, 6};
	double z[6] = {1, -2, 1, 4, 1, 6};
	double w[6] = {1, -2, 1, 4, 1, 6};
	double B[9] = {14, 23, 16, 0, 38, 26, 0, 0, 35};
	double C[9] = {14, 23, 16, 0, 38, 26, 0, 0, 35};
	double det_a;
	int p[3] = {0};
	int n = 3;
	int nrhs = 1;

	printf( "A=\n" );
	dense_print_matrix( n, n, A, type );
	det_a = dense_eval_det( n, A, FACTOR_LU_RIGHT_LOOKING, type );
	printf( "\n|A|=%.10le\n", det_a );
	printf( "\nx=\n" );
	dense_print_vector( n, x, type );
	printf( "\ny=\n" );
	dense_print_vector( 2 * n, y, type );
	printf( "\nB=\n" );
	dense_print_matrix_trig( n, B, TRIG_LOWER, type );
	printf( "\nz=\n" );
	dense_print_vector( 2 * n, z, type );

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
	det_a = dense_eval_factor_det( n, A, p, FACTOR_LU_RIGHT_LOOKING, type );
	printf( "\n|A|=%.10le\n", det_a );

	dense_solve( n, nrhs, A, x, p, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, type );
	printf( "\nsolve x=\n" );
	dense_print_vector( n, x, type );

	nrhs = 2;
	dense_solve( n, nrhs, A, y, p, FACTOR_LU_RIGHT_LOOKING, TRANS_NONE, type );
	printf( "\nsolve y1=\n" );
	dense_print_vector( n, y, type );
	printf( "\nsolve y2=\n" );
	dense_print_vector( n, y + n, type );

	if ( !dense_lu_factor( n, B, NULL, FACTOR_LU_CHOLESKY, type ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}
	printf( "\nLU=\n" );
	dense_print_matrix_trig( n, B, TRIG_LOWER, type );
	nrhs = 2;
	dense_solve( n, nrhs, B, z, NULL, FACTOR_LU_CHOLESKY, TRANS_NONE, type );
	printf( "\nsolve z1=\n" );
	dense_print_vector( n, z, type );
	printf( "\nsolve z2=\n" );
	dense_print_vector( n, z + n, type );

	if ( !dense_lu_factor( n, C, p, FACTOR_LU_BUNCH_KAUFMAN, type ) )
	{
		fprintf( stderr, "[Error] LU factorization fail\n" );
		abort();
	}
	printf( "\nLU=\n" );
	dense_print_matrix_trig( n, C, TRIG_LOWER, type );
	printf( "\nperm=\n" );
	dense_print_vector_i( n, p );
	printf( "\nP=\n" );
	dense_print_matrix_perm( n, p );
	nrhs = 2;
	dense_solve( n, nrhs, C, w, p, FACTOR_LU_BUNCH_KAUFMAN, TRANS_NONE, type );
	printf( "\nsolve w1=\n" );
	dense_print_vector( n, w, type );
	printf( "\nsolve w2=\n" );
	dense_print_vector( n, w + n, type );

	return EXIT_SUCCESS;
}



