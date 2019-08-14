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

void FUNC_WITH_POSTFIX( interpolation_lagrange, POSTFIX ) (
	long n, // Pn(xp) = yp
	long ldy, // length of yp
	double *x, 
	double **y,
	double xp,
	double *yp
);

#endif

