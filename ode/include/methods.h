#ifndef METHODS_H
#define METHODS_H

// Adams-Moulton
// multistep implicit method
// yn = yn_1 + integral(P(x)) from xn_1 to xn
// P(x) is interpolation polynomil of f(x,y(x)) from xn, xn_1, ...
double adams_moulton ( int order, double xn_1, double yn_1, double hn_1, double *ylist );

// Adams-Bashforth
// multistep explicit method
// yn = yn_1 + integral(P(x)) from xn_1 to xn
// P(x) is interpolation polynomil of f(x,y(x)) from xn_1, xn_2, ...
double adams_bashforth ( int order, double xn_1, double yn_1, double hn_1, double *ylist );

// Backward Differential Formula (Gear)
// multistep implicit method
// yn = yn_1 + integral(P(x)) from xn_1 to xn
// P(x) is interpolation polynomil of f(x,y(x)) from xn, xn_1, ...
double bdf ( int order, double xn_1, double yn_1, double hn_1, double *ylist );

// other...
extern double diff  ( double yn, double xn );
extern double diff2 ( double yn, double xn );
extern double diff_exact ( double yn, double xn );
extern double exact ( double xn );

extern double g_initial_solution;

#endif
