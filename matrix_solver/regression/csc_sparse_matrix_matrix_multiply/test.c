#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	g_matrix_print_format = MATRIX_PRINT_FORMAT_MATLAB;

	sparse_csc_t A;
	A.nz = 9;
	A.m = 4;
	A.n = 3;
	A.Ap = (sparse_int *) calloc ( A.n + 1, sizeof(sparse_int) );
	A.Ai = (sparse_int *) calloc ( A.nz, sizeof(sparse_int) );
	A.Ax = (sparse_float *) calloc ( A.nz, sizeof(sparse_float) );
	A.xtype = REAL_NUMBER;

	A.Ap[0] = 0;
	A.Ai[0] = 0;
	A.Ax[0] = 4.1;
	A.Ai[1] = 1;
	A.Ax[1] = 3.2;
	A.Ai[2] = 2;
	A.Ax[2] = 2.3;
	A.Ai[3] = 3;
	A.Ax[3] = 1.4;
	A.Ap[1] = 4;

	A.Ai[4] = 1;
	A.Ax[4] = 0.2;
	A.Ai[5] = 3;
	A.Ax[5] = 0.1;
	A.Ap[2] = 6;

	A.Ai[6] = 0;
	A.Ax[6] = -1;
	A.Ai[7] = 2;
	A.Ax[7] = -2;
	A.Ai[8] = 3;
	A.Ax[8] = -3;
	A.Ap[3] = 9;

	printf( "A = " );
	sparse_print_csc_full_matrix ( &A );


	sparse_triplet_t B;
	B.nz = 7;
	B.m = 3;
	B.n = 4;
	B.elements = (sparse_element_t *) calloc ( B.nz, sizeof(sparse_element_t) ) ;
	B.elements[0].row = 2;
	B.elements[0].col = 2;
	B.elements[0].x = 3.1;
	B.elements[1].row = 0;
	B.elements[1].col = 0;
	B.elements[1].x = 1.2;
	B.elements[2].row = 0;
	B.elements[2].col = 1;
	B.elements[2].x = 100;
	B.elements[3].row = 0;
	B.elements[3].col = 3;
	B.elements[3].x = 200;
	B.elements[4].row = 1;
	B.elements[4].col = 2;
	B.elements[4].x = 300;
	B.elements[5].row = 2;
	B.elements[5].col = 1;
	B.elements[5].x = -20;
	B.elements[6].row = 2;
	B.elements[6].col = 3;
	B.elements[6].x = -10;
	printf( "B = " );

	sparse_csc_t *B_csc = sparse_convert_triplet_to_CSC( &B );
	sparse_print_csc_full_matrix ( B_csc );
	sparse_csc_t *C = sparse_matrix_matrix_multiply ( &A, B_csc );

	printf( "C = " );
	sparse_print_csc_full_matrix ( C );

	return EXIT_SUCCESS;
}

