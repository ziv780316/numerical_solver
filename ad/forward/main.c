#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <complex.h>

double f3_r_diff ( double x )
{
	return 1.0 / x;
}

double f3_r ( double x )
{
	return log(x);
}

complex f3_z ( complex z )
{
	return clog(z);
}

double f2_r_diff ( double x )
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

double f1_r_diff ( double x )
{
	return 2.0 * x;
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

// Error = O(h^2) and does no subject to substract error when h is very small
// f(x + ih) ~= f(x) + ih*f^1(x) - h^2*f^2(x) + ih^3*f^3(x) -...
// f^1(x) = Im(f(x + ih))/h, when h -> 0 
double complex_ad ( double x, double h, complex (*f) (complex) )
{
	complex zin  = x + (I * h);
	complex zout = f( zin );
	return cimag(zout) / h;
}

int main ( int argc __attribute__((unused)), char **argv __attribute__((unused)) )
{
	double x = 1e-3;
	double h = 1.0;
	double (*f_r) (double) = f3_r;
	complex (*f_z) (complex) = f3_z;
	double exact = f3_r_diff(x);
	double fd;
	double cd;
	double ad;
	double error_fd;
	double error_cd;
	double error_ad;

	printf( "# DBL_EPSILON=%.10e\n", DBL_EPSILON );
	printf( "# exact=%.10e\n", exact );
	printf( "# h fd cd ad\n" );

	for ( int i = 1; i <= 20; ++i )
	{
		h *= 0.1;

		fd = forward_difference( x, h, f_r );
		cd = central_difference( x, h, f_r );
		ad = complex_ad( x, h, f_z );

		// finite difference with be 0 when h is very small (round-off error due to substraction)
		error_fd = fabs((exact - fd) / exact);
		error_cd = fabs((exact - cd) / exact);
		error_ad = fabs((exact - ad) / exact); 

		printf( "%.10e %.10e %.10e %.10e\n", h, error_fd, error_cd, error_ad );
	}

	return EXIT_SUCCESS;
}

