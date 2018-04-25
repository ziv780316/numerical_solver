#ifndef METHODS_H
#define METHODS_H

// Adams-Moulton
// multistep implicit method
// yn = yn_1 + integral(P(t)) from tn_1 to tn
// P(t) is interpolation polynomil of f(t,y(t)) from tn, tn_1, ...
double adams_moulton ( int order, double tn_1, double yn_1, double hn_1, double *ylist );

// Adams-Bashforth
// multistep explicit method
// yn = yn_1 + integral(P(t)) from tn_1 to tn
// P(t) is interpolation polynomil of f(t,y(t)) from tn_1, tn_2, ...
double adams_bashforth ( int order, double tn_1, double yn_1, double hn_1, double *ylist );

// Backward Differential Formula (Gear)
// multistep implicit method
// yn = yn_1 + integral(P(t)) from tn_1 to tn
// P(t) is interpolation polynomil of f(t,y(t)) from tn, tn_1, ...
double bdf ( int order, double tn_1, double yn_1, double hn_1, double *ylist );

// other...
extern double diff  ( double yn, double tn );
extern double exact ( double tn );
double jacobian ( double y, double t );

extern double g_initial_solution;
extern int g_total_iteration;

#endif
