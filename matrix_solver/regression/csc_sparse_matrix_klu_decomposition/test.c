#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "matrix_solver.h"

int main ( int argc, char **argv )
{
	g_matrix_print_format = MATRIX_PRINT_FORMAT_MATLAB;
	g_debug_sparse_lu_decomposition = true;

	sparse_triplet_t A;
	A.nz = 13;
	A.m = 4;
	A.n = 4;
	A.elements = (sparse_element_t *) calloc ( A.nz, sizeof(sparse_element_t) ) ;
	A.elements[0].row = 0;
	A.elements[0].col = 0;
	A.elements[0].x = 1;
	A.elements[1].row = 1;
	A.elements[1].col = 0;
	A.elements[1].x = 2;
	A.elements[2].row = 2;
	A.elements[2].col = 0;
	A.elements[2].x = 3;
	A.elements[3].row = 3;
	A.elements[3].col = 0;
	A.elements[3].x = 4;
	A.elements[4].row = 0;
	A.elements[4].col = 1;
	A.elements[4].x = 2;
	A.elements[5].row = 2;
	A.elements[5].col = 1;
	A.elements[5].x = 5;
	A.elements[6].row = 3;
	A.elements[6].col = 1;
	A.elements[6].x = -10;
	A.elements[7].row = 0;
	A.elements[7].col = 2;
	A.elements[7].x = 1;
	A.elements[8].row = 1;
	A.elements[8].col = 2;
	A.elements[8].x = 6;
	A.elements[9].row = 3;
	A.elements[9].col = 2;
	A.elements[9].x = -7;
	A.elements[10].row = 0;
	A.elements[10].col = 3;
	A.elements[10].x = 2.2;
	A.elements[11].row = 1;
	A.elements[11].col = 3;
	A.elements[11].x = 5.3;
	A.elements[12].row = 2;
	A.elements[12].col = 3;
	A.elements[12].x = -1.8;

	long p_user[4];
	p_user[0] = 1;
	p_user[1] = 0;
	p_user[2] = 2;
	p_user[3] = 3;


	printf( "A = " );
	sparse_csc_t *A_csc = sparse_convert_triplet_to_CSC( &A );
	sparse_print_csc_full_matrix ( A_csc );

	sparse_csc_t *P, *L, *U, *Q, *R;
	sparse_float *r;
	sparse_int *p, *pinv, *q;
	if ( !sparse_matrix_lu_decomposition ( A_csc, SPARSE_LU_DECOMPOSITION_KLU, &L, &U, &P, &Q, &R, &p, &pinv, &q, &r, p_user, NULL, 0, 1 ) )
	{
		fprintf( stderr, "[Error] KLU decomposition fail\n" );
		exit(1);
	}

	// P*R*A*Q = L*U

	printf( "R = " );
	sparse_print_csc_full_matrix ( R );

	printf( "Rinv = " );
	sparse_matrix_diagonal_inverse ( R );
	sparse_print_csc_full_matrix ( R );

	printf( "1/Rs = " );
	dense_print_vector ( A.n, r, REAL_NUMBER );

	printf( "P = " );
	sparse_print_csc_full_matrix ( P );

	printf( "p = [...\n" );
	for ( int i = 0; i < 4; ++i )
	{
		printf( "%ld\n", p[i] );
	}
	printf( "];\n" );

	printf( "pinv = [...\n" );
	for ( int i = 0; i < 4; ++i )
	{
		printf( "%ld\n", pinv[i] );
	}
	printf( "];\n" );

	printf( "L = " );
	sparse_print_csc_full_matrix ( L );

	printf( "U = " );
	sparse_print_csc_full_matrix ( U );

	printf( "Q = " );
	sparse_print_csc_full_matrix ( Q );

	printf( "Pinv = " );
	sparse_matrix_transpose ( P );
	sparse_print_csc_full_matrix ( P );

	printf( "Qinv = " );
	sparse_matrix_transpose ( Q );
	sparse_print_csc_full_matrix ( Q );

	sparse_csc_t *L_U;
	printf( "L_U = " );
	L_U = sparse_matrix_matrix_multiply ( L, U );
	sparse_print_csc_full_matrix ( L_U );

	sparse_csc_t *inv_P_L_U;
	printf( "inv_P_L_U = " );
	inv_P_L_U = sparse_matrix_matrix_multiply ( P, L_U );
	sparse_print_csc_full_matrix ( inv_P_L_U );

	sparse_csc_t *inv_R_inv_P_L_U_inv_Q;
	inv_R_inv_P_L_U_inv_Q = sparse_matrix_matrix_multiply ( R, P );
	printf( "inv_R_inv_P = " );
	sparse_print_csc_full_matrix ( inv_R_inv_P_L_U_inv_Q );
	inv_R_inv_P_L_U_inv_Q = sparse_matrix_matrix_multiply ( inv_R_inv_P_L_U_inv_Q, L );
	inv_R_inv_P_L_U_inv_Q = sparse_matrix_matrix_multiply ( inv_R_inv_P_L_U_inv_Q, U );
	inv_R_inv_P_L_U_inv_Q = sparse_matrix_matrix_multiply ( inv_R_inv_P_L_U_inv_Q, Q );
	printf( "inv_R_inv_P_L_U_inv_Q = " );
	sparse_print_csc_full_matrix ( inv_R_inv_P_L_U_inv_Q );

	return EXIT_SUCCESS;
}

