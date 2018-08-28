#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "newton.h"
#include "matrix_solver.h"


static void check_user_define_args ( double *x0, 
				     void (load_f) (double *x, double*f),
				     void (load_jacobian) (double *x, double*J) )
{
	if ( !x0 )
	{
		printf( "[Error] x0 is undefined\n" );
		abort();
	}
	if ( !load_f )
	{
		printf( "[Error] cannot find F(x)\n" );
		abort();
	}
	if ( !load_jacobian )
	{
		printf( "[Warning] cannot find F'(x), use approximate derivative by finite difference\n" );
	}
}

static void newton_initialize ( int n, double *x, double *x0, bool random_initial )
{
	memcpy( x, x0, sizeof(double) * n );	
}

// J = J + (df - J*dx)*dxT / |dx|^2
static void broyden_update ( int n, double *J, double *df, double *dx, bool debug )
{
	double dx_square;
	double *work = (double *) malloc ( sizeof(double) * n );
	memcpy( work, df, sizeof(double) * n );

	// |dx|^2
	dense_vector_inner_product ( n, dx, dx, &dx_square );

	// work = df - (J * dx)
	dense_matrix_vector_multiply ( n, n, -1.0, J, dx, 1.0, work, false );

	// J = J + (x . yT) / |dx|^2
	dense_maxtrix_rank_1_update ( n, J, 1.0/dx_square, work, dx );

	free( work );
}

// J = J + (dx - J*df)*dxT*J / (dxT*J*df)
static void broyden_update_sherman_morrison ( int n, double *J, double *df, double *dx, bool debug )
{
	double dx_square;
	double *work1 = (double *) malloc ( sizeof(double) * n );
	double *work2 = (double *) malloc ( sizeof(double) * n );

	// work = dx - (J * df)
	memcpy( work1, dx, sizeof(double) * n );
	dense_matrix_vector_multiply ( n, n, -1.0, J, df, 1.0, work1, false );

	// dxT*J = JT*dx
	dense_matrix_vector_multiply ( n, n, 1.0, J, dx, 0.0, work2, true );

	// (dxT*J) * df
	dense_vector_inner_product ( n, work2, df, &dx_square );

	// J = J + (x . yT) / |dx|^2
	dense_maxtrix_rank_1_update ( n, J, 1.0/dx_square, work1, work2 );

	free( work1 );
	free( work2 );
}

bool newton_solve ( newton_iterative_type iterative_type, 
		    newton_modified_type modified_type,
		    newton_derivative_type diff_type,
		    int n,
		    double *x0,
		    double *x_result,
		    double *f_result,
		    void (load_f) (double *x, double*f),
		    void (load_jacobian) (double *x, double*J),
		    int maxiter,
		    int miniter,
		    int *total_iter,
		    double rtol,
		    double atol,
		    double residual_tol,
		    bool random_initial,
		    bool debug )
{
	// check load_f exist or not
	check_user_define_args( x0, load_f, load_jacobian );

	// initialize first solution
	int J_size = n * n;
	int *perm = (int *) malloc ( sizeof(int) * n ); // use in lapack LU factorization
	double *x = (double *) malloc ( sizeof(double) * n );
	double *dx = (double *) malloc ( sizeof(double) * n );
	double *x_candidate = NULL;
	double *f = (double *) malloc ( sizeof(double) * n );
	double *df = NULL;
	double *f_candidate = NULL;
	double *f_delta_forward = NULL;
	double *f_delta_backward = NULL;
	double *rhs = (double *) malloc ( sizeof(double) * n );
	double *J = (double *) malloc ( sizeof(double) * J_size );
	double *J_old = NULL;

	if ( diff_type != NEWTON_DIFF_JACOBIAN )
	{
		f_delta_forward  = (double *) malloc ( sizeof(double) * n );
		if ( diff_type == NEWTON_DIFF_CENTRAL )
		{
			f_delta_backward = (double *) malloc ( sizeof(double) * n );
		}
	}
	if ( NEWTON_BROYDEN == iterative_type )
	{
		df = (double *) malloc ( sizeof(double) * n );
		J_old = (double *) malloc ( sizeof(double) * J_size );
	}
	else if ( (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		  (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
	{
		df = (double *) malloc ( sizeof(double) * n );
	}
	if ( modified_type != MODIFIED_NONE )
	{
		x_candidate = (double *) malloc ( sizeof(double) * n );
		f_candidate = (double *) malloc ( sizeof(double) * n );
	}

	newton_initialize( n, x, x0, random_initial );

	// iterative procedure
	int iter = 1;
	bool converge = false;
	double delta_ratio = 1e-9;
	double delta;
	double delta_inv;
	double x_tmp;
	while ( !converge )
	{
		if ( (iter > maxiter) && (-1 != maxiter) && (iter >= miniter) )
		{
			break;
		}

		// load RHS
		if ( (NEWTON_BROYDEN == iterative_type) ||
		     (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		     (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
		{
			memcpy( df, f, sizeof(double) * n );	
		}
		load_f( x, f );
		memcpy( rhs, f, sizeof(double) * n );	
		if ( (NEWTON_BROYDEN == iterative_type) ||
		     (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		     (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
		{
			for ( int i = 0; i < n; ++i )
			{
				df[i] = f[i] - df[i];
			}
		}
		if ( (1 == iter) && debug )
		{
			printf( "------- initial -------\n" );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d]=%.10e f=%.10e\n", i, x[i], f[i] );
			}
		}

		// construct jacobian matrix
		if ( (NEWTON_NORMAL == iterative_type) || 
		     ((NEWTON_CHORD == iterative_type) && (1 == iter)) ||
		     ((NEWTON_BROYDEN == iterative_type) && (1 == iter)) ||
		     ((NEWTON_BROYDEN_INVERTED == iterative_type) && (1 == iter)) ||
		     ((NEWTON_BROYDEN_INVERTED_BAD == iterative_type) && (1 == iter)) 
		   )
		{
			if ( (NEWTON_DIFF_JACOBIAN  == diff_type) && load_jacobian )
			{
				// use user pre-define jacobian 
				load_jacobian( x, J );
			}
			else
			{
				for ( int i = 0; i < n; ++i )
				{
					if ( fabs(x[i]) > DBL_EPSILON )
					{
						if ( NEWTON_DIFF_FORWARD == diff_type )
						{
							// use forward difference for better speed
							delta = x[i] * delta_ratio;
							delta_inv = 1.0 / delta;
							x_tmp = x[i];
							x[i] += delta;
							load_f( x, f_delta_forward );
							x[i] = x_tmp;
							for ( int k = 0; k < n; ++k )
							{
								*(J + n*i + k) = (f_delta_forward[k] - f[k]) * delta_inv;
							}
						}
						else if ( NEWTON_DIFF_CENTRAL == diff_type )
						{
							// use central difference for accurate derivative approximation O(h^2)
							delta = x[i] * delta_ratio;
							delta_inv = 1.0 / (2.0 * delta);
							x_tmp = x[i];
							x[i] += delta;
							load_f( x, f_delta_forward );
							x[i] = x_tmp;
							x[i] -= delta;
							load_f( x, f_delta_backward );
							x[i] = x_tmp;
							for ( int k = 0; k < n; ++k )
							{
								*(J + n*i + k) = (f_delta_forward[k] - f_delta_backward[k]) * delta_inv;
							}
						}
					}
					else
					{
						if ( debug )
						{
							printf( "iter=%d x[%d]=%.10e near zero, skip finite difference\n", iter, i, x[i] );
						}
					}
				}
			}

			if ( (NEWTON_BROYDEN == iterative_type) && (1 == iter) )
			{
				// save old J before factorization
				memcpy( J_old, J, sizeof(double) * J_size );	

			}
		}
		else if ( ((NEWTON_BROYDEN == iterative_type) && (iter > 1)) )
		{
			// Broyden update J
			broyden_update( n, J_old, df, dx, debug );
			memcpy( J, J_old, sizeof(double) * J_size );	
		}

		if ( debug )
		{
			if ( (NEWTON_BROYDEN_INVERTED != iterative_type) &&  
			     (NEWTON_BROYDEN_INVERTED_BAD != iterative_type) )
			{
				printf( "J = \n" );
				dense_print_matrix ( n, n, J );
			}
		}

		// matrix factorization A = PLU
		int matrix_ok;
		if ( !((NEWTON_CHORD == iterative_type) && (iter > 1)) &&
		     !((NEWTON_BROYDEN_INVERTED == iterative_type) && (iter > 1)) &&
		     !((NEWTON_BROYDEN_INVERTED_BAD == iterative_type) && (iter > 1)) )
		{
			matrix_ok = dense_lu_factor ( n, J, perm );
			if ( !matrix_ok )
			{
				fprintf( stderr, "[Error] LU factorization fail\n" );
				abort();
			}
		}

		// solve J(dx) = -F
		for ( int i = 0; i < n; ++i )
		{
			rhs[i] *= -1.0;
		}
		if ( (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		     (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
		{
			if ( 1 == iter )
			{
				matrix_ok = dense_matrix_inverse ( n, J, perm );
				if ( !matrix_ok )
				{
					fprintf( stderr, "[Error] inverse jacobian matrix fail\n" );
					abort();
				}
			}
			else 
			{
				
				if ( NEWTON_BROYDEN_INVERTED == iterative_type )
				{
					// Sherman-Morrison Broyden update, the same of bad Broyden but with better numerical stability
					broyden_update_sherman_morrison( n, J, df, dx, debug );
				}
				else if ( NEWTON_BROYDEN_INVERTED_BAD == iterative_type )
				{
					// bad Broyden update, more efficient but not numerical stable 
					broyden_update( n, J, dx, df, debug );
				}
			}

			if ( debug )
			{
				dense_lu_factor ( n, J, perm );
				dense_matrix_inverse ( n, J, perm );
				printf( "J = \n" );
				dense_print_matrix ( n, n, J );

				dense_lu_factor ( n, J, perm );
				dense_matrix_inverse ( n, J, perm );
				printf( "J^-1 = \n" );
				dense_print_matrix ( n, n, J );
			}

			// dx = J^-1 * -f
			dense_matrix_vector_multiply ( n, n, 1.0, J, rhs, 0.0, dx, false );
		}
		else
		{
			matrix_ok = dense_solve ( n, J, rhs, perm, false );
			if ( !matrix_ok )
			{
				fprintf( stderr, "[Error] LU solve fail\n" );
				abort();
			}
			memcpy( dx, rhs, sizeof(double) * n );
		}

		// show solve results 
		if ( debug )
		{
			printf( "------- iter %d -------\n", iter );
			for ( int i = 0; i < n; ++i )
			{
				printf( "x[%d] new=%.10e old=%.10e dx=%.10e f=%.10e\n", i, x[i] + dx[i], x[i], dx[i], f[i] );
			}
		}

		// check stop criteria
		double diff;
		double tol;
		converge = true;
		for ( int i = 0; i < n; ++i )
		{
			// check iteration converge
			diff = abs( dx[i] );
			tol = abs(x[i]) * rtol + atol;
			if ( diff > tol )
			{
				if ( debug )
				{
					printf( "iter=%d x[%d]=%.10e x_new[%d]=%.10e does not converged, diff=%.10e, tol=%.10e\n", iter, i, x[i], i, x[i] + dx[i], diff, tol );
				}
				converge = false;
				break;
			}

			// check residual converge (maximum norm)
			if ( fabs(f[i]) > residual_tol )
			{
				if ( debug )
				{
					printf( "iter=%d f[%d]=%.10e does not converged residual_tol=%.10e\n", iter, i, f[i], residual_tol );
				}
				converge = false;
				break;
			}
		}

		// modified newton for better performance or prevent too large step
		if ( MODIFIED_LINE_SEARCH == modified_type )
		{
		}

		// next iteration
		for ( int i = 0; i < n; ++i )
		{
			x[i] += dx[i];
		}

		// for benchmark performance
		if ( iter < miniter )
		{
			converge = false;
		}

		++iter;
		++(*total_iter);
	}

	// store final results
	memcpy( x_result, x, sizeof(double) * n );	
	memcpy( f_result, f, sizeof(double) * n );	

	// show newton results
	if ( debug )
	{
		printf( "\n========== Newton Converge %s in %d Iteration ==========\n", (converge ? "Success" : "Fail"), iter - 1 );
		for ( int i = 0; i < n; ++i )
		{
			printf( "x[%d]=%.10e  f=%.10e\n", i, x[i], f[i] );
		}
	}

	// release memory
	free( perm );
	free( x );
	free( dx );
	free( f );
	free( rhs );
	free( J );
	if ( diff_type != NEWTON_DIFF_JACOBIAN )
	{
		free( f_delta_forward  );
		if ( diff_type == NEWTON_DIFF_CENTRAL )
		{
			free( f_delta_backward );
		}
	}
	if ( NEWTON_BROYDEN == iterative_type )
	{
		free( df );
		free( J_old );
	}
	else if ( (NEWTON_BROYDEN_INVERTED == iterative_type) || 
		  (NEWTON_BROYDEN_INVERTED_BAD == iterative_type) )
	{
		free( df );
	}
	if ( modified_type != MODIFIED_NONE )
	{
		free( x_candidate );
		free( f_candidate );
	}

	return converge;
}
