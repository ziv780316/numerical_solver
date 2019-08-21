#ifndef NUMERICAL_UTIL_H
#define NUMERICAL_UTIL_H

#if defined PRECISION_32_BIT
#define POSTFIX 32
#define FLOAT_TYPE float

#elif defined PRECISION_64_BIT
#define POSTFIX 64
#define FLOAT_TYPE double

#elif defined PRECISION_80_BIT
#define POSTFIX 80
#define FLOAT_TYPE long double

#elif defined PRECISION_128_BIT
#define POSTFIX 128
#define FLOAT_TYPE __float128

#else
#define POSTFIX 64
#define FLOAT_TYPE double

#endif

#define __FUNC_WITH_POSTFIX(func, postfix) func ## _ ## postfix
#define FUNC_WITH_POSTFIX(func, postfix) __FUNC_WITH_POSTFIX(func, postfix)

// efficient for interpolation multiple Y
// Pn(xp) = l0*y0 + l1*y1 + l2*y2 + ... + ln*yn
void FUNC_WITH_POSTFIX( interpolation_lagrange, POSTFIX ) (
	long n, // Pn(xp) = yp
	long ldy, // length of yp
	double *x, 
	double **y,
	double xp,
	double *yp
);

// backward difference y0 = y[0], yn = y[n]
// ddN = dd[n], dd1 = dd[1], y0 = dd[0]
// ddN = [yn, yn-1, ..., y0] = (y⁽ⁿ⁾(ζ) / N!) where xn ≤ ζ ≤ x0
void FUNC_WITH_POSTFIX( divide_difference, POSTFIX ) (
	long n, // backward difference [yn, yn-1, ..., y0], i.e. (y1 - y0)/(x1 - x0) and x0 > x1
	double *x, // x0 is x[0], xn is x[n]
	double *y, // backward difference, y0 is y[0], yn is y[n]
	double *dd 
);

// efficient for interpolation multiple X
// ddN = dd[n], dd1 = dd[1], y0 = dd[0]
// Pn(xp) = dd[0] + dd[1]*(x-x0) + dd[2]*(x-x0)*(x-x1) + ... + dd[n]*(x-x0)*...(x-xn-1)
void FUNC_WITH_POSTFIX( interpolation_newton, POSTFIX ) (
	long n, 
	int ldx,
	double *x, 
	double *xp, 
	double *yp, 
	double *dd 
);

#endif

