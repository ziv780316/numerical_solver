#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <complex.h>

double f3_r_diff ( double x, int order )
{
	double diff = 0.0;

	switch ( order )
	{
		case 1:
			diff = 1.0 / x;
		break;

		case 2:
			diff = -1.0 / (x*x);
		break;
	}

	return diff;
}

double f3_r ( double x )
{
	return log(x);
}

complex f3_z ( complex z )
{
	return clog(z);
}

double f2_r_diff ( double x, int order __attribute__((unused)) )
{
	return exp(x);
}

double f2_r ( double x )
{
	return exp(x);
}

complex f2_z ( complex z )
{
	return cexp(z);
}

double f1_r_diff ( double x, int order )
{
	double diff = 0.0;

	switch ( order )
	{
		case 1:
			diff = 2.0 * x;
		break;

		case 2:
			diff = 2.0;
		break;
	}

	return diff;
}

double f1_r ( double x )
{
	return pow(x, 2.0);
}

complex f1_z ( complex z )
{
	return cpow(z, 2.0);
}

// Error = O(h)
double forward_difference ( double x, double h, double (*f) (double) )
{
	return (f(x + h) - f(x)) / h;
}

// Error = O(h^2) but subject to substract error when h is very small
double central_difference ( double x, double h, double (*f) (double) )
{
	return (f(x + h) - f(x - h)) / (2.0 * h);
}

// first order derivative error = O(h^2) and does no subject to substract error when h is very small
// f(x + ih) ~= f(x) + ih*f^1(x) - h^2*f^2(x)/2 + ih^3*f^3(x)/6 -...
// f^1(x) = Im(f(x + ih))/h, when h -> 0 
// f^2(x) = (f(x) - Re(f(x + ih)))/(0.5*(h^2)), when h -> 0 
double complex_ad ( double x, double h, complex (*f) (complex), int order )
{
	double diff = 0.0;
	complex zin;
	complex zout;

	switch ( order )
	{
		case 1:
			// O(h^2) without round-off error of substraction
			zin  = x + (I * h);
			zout = f( zin );
			diff = cimag(zout) / h;
			break;

		case 2:
			// O(h^4) but subject to round-off error
			zin  = x + (I * h);
			zout = f( zin );
			diff = (f(x) - creal(zout)) / (0.5*h*h);
			break;

		default:
			fprintf( stderr, "[Error] cannot support order %d derivative with finite difference\n", order );
			break;
	}

	return diff;
}

// finite difference for high order derivative
double finite_difference ( double x, double h, double (*f) (double), int order )
{
	double diff = 0.0;

	switch ( order )
	{
		case 1:
			// O(h^2) but subject to round-off error
			diff = central_difference( x, h, f ); 
			break;

		case 2:
			// O(h^4) but subject to round-off error
			diff = (-f(x + 2*h) + 16*f(x + h) - 30*f(x) + 16*f(x - h) - f(x - 2*h)) / (12*h*h); 
			break;

		default:
			fprintf( stderr, "[Error] cannot support order %d derivative with finite difference\n", order );
			break;
	}

	return diff;
}

int main ( int argc __attribute__((unused)), char **argv __attribute__((unused)) )
{
	double x = 1e-3;
	double h = 1.0;
	double (*f_r) (double) = f3_r;
	complex (*f_z) (complex) = f3_z;
	double exact = f3_r_diff( x, 1 );
	double exact2 = f3_r_diff( x, 2 );
	double fd;
	double fd2;
	double cd;
	double ad;
	double ad2;
	double error_fd;
	double error_fd2;
	double error_cd;
	double error_ad;
	double error_ad2;

	printf( "# DBL_EPSILON=%.10e\n", DBL_EPSILON );
	printf( "# exact=%.10e\n", exact );
	printf( "# h fd cd ad fd2 ad2\n" );

	for ( int i = 1; i <= 20; ++i )
	{
		h *= 0.1;

		fd = forward_difference( x, h, f_r );
		cd = central_difference( x, h, f_r );
		ad = complex_ad( x, h, f_z, 1 );

		fd2 = finite_difference( x, h, f_r, 2 );
		ad2 = complex_ad( x, h, f_z, 2 );

		// finite difference with be 0 when h is very small (round-off error due to substraction)
		error_fd = fabs((exact - fd) / exact);
		error_cd = fabs((exact - cd) / exact);
		error_ad = fabs((exact - ad) / exact); 
		error_fd2 = fabs((exact2 - fd2) / exact2);
		error_ad2 = fabs((exact2 - ad2) / exact2); 

		printf( "%.10e %.10e %.10e %.10e %.10e %.10e\n", h, error_fd, error_cd, error_ad, error_fd2, error_ad2 );
	}

	return EXIT_SUCCESS;
}

