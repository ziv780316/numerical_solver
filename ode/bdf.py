#!/usr/bin/env python3.8
import math
import copy
import numpy as np

# ---------------------------------------------------------------------------------
# Numerical Routine
# ---------------------------------------------------------------------------------

# numerical discretization method
def discretize_bdf ( order, h_list ):
  '''
  BDF - Backward Differentiation Formula
  yₙ₊₁ = α₀yₙ + α₁yₙ₋₁ + α₂yₙ₋₂ + ... αₖyₙ₋ₖ+ β₀y'ₙ₊₁ 
  order = k + 1
  '''
  debug = True
  n = order + 1
  A = np.zeros( (n,n) )
  b = np.zeros( (n,1) )
  for i in range(0, n-1):
    A[0, i] = 1
  for i in range(0, n):
    b[i] = (h_list[0]**i) / math.factorial(i)
  for i in range(1, n):
    A[i,n-1] = (h_list[0]**i) / math.factorial(i-1)
  for i in range(1, n):
    h_sum = 0
    for j in range(1, n-1):
      h_sum += h_list[j]
      A[i, j] = ((-h_sum)**i) / math.factorial(i)

  if debug:
    print( 'A=\n', A );
    print( 'b=\n', b );

  c = np.linalg.solve(A, b)

  if debug:
    for i in range(0, order):
      print( 'h{id} = {h}'.format(id=i, h=h_list[i]) );
    for i in range(1, n):
      print( '\u03b1{id} = {c}'.format(id=i-1, c=c[i-1]) );
    print( '\u03b20 = {c}'.format(c=c[n-1,0]) );

  return c

def discretize_bdf_diff ( order, h_list ):
  debug = True

  n = order + 1
  c = discretize_bdf( order, h_list )
  beta = c[order]
  c_diff = np.zeros( (n) )
  c_diff[0] = 1/beta
  for i in range(1, n):
    c_diff[i] = -c[i-1]/beta

  if debug:
    for i in range(0, n):
      print( '\u03b1_diff{id} = {c}'.format(id=i, c=c_diff[i]) );

  return c_diff

# interpolation polynomial
def interp_poly ( order, t_list, t ):
  '''
  Lagrange Interpolation
  '''
  debug = True
  n = order + 1
  l = np.ones( (n) )
  for i in range(0, n):
    for j in range(0, n):
      if j == i:
        continue 
      l[i] *= (t - t_list[j])/(t_list[i]-t_list[j])

  if debug:
    for i in range(0, n):
      print( 'L{id} = {c}'.format(id=i, c=l[i]) );

  return l

def interp_diff_poly ( order, t_list, t ):
  '''
  Lagrange Differentiation Interpolation
  '''
  debug = True
  n = order + 1
  ld = np.zeros( (n) )
  for i in range(0, n):

    denom = 1
    for j in range(0, n):
      if j == i:
        continue 
      denom *= (t_list[i]-t_list[j])

    numer = 0
    for j in range(0, n):
      if j == i:
        continue 
      tmp = 1
      for k in range(0, n):
        if k == i or k == j:
          continue 
        tmp *= (t - t_list[k])
      numer += tmp
    ld[i] = numer / denom

  if debug:
    for i in range(0, n):
      print( 'Ld{id} = {c}'.format(id=i, c=ld[i]) );

  return ld

# interpolation by polynomial
def interp ( l, x ):
  '''
  Lagrange Interpolation
  '''
  xp = np.dot(l, x)
  return xp

def get_t_list (t, h_list):
  n = max( h_list.shape ) + 1
  t_list = np.zeros( (n) )
  t_list[0] = t
  for i in range(1, n):
    t_list[i] = t_list[i-1] - h_list[i-1]

  return t_list

  
# NR iteration
def do_nr ( t, dt, x, xd ):
  rtol = 1e-3
  atol = 1e-9
  ftol = 1e-9

  iter_max = 10
  iterno = 1
  while iterno < iter_max:
    xd = discretize_bdf( order, h_list ) 
    f = test_f1( t, dt, x, xd )
    iterno += 1


# ---------------------------------------------------------------------------------
# User Routine
# ---------------------------------------------------------------------------------

# test problem
def test_f1 ( t, dt, x, xd ):
  lamda = 1
  f = (x + dt*lamda*xd) 
  return f


# test function 1
def exact_x_problem_1( t ):
  return t**2 + t

def exact_xd_problem_1( t ):
  return 2*t + 1

def test_problem_1 ():
  '''
  x(t) = t^2 + t
  x'(t) = 2*t + 1
  '''
  order = 3
  n = order + 1 # time point number
  h_list = np.zeros( (order) )
  h_list[0] = 1
  h_list[1] = 2
  h_list[2] = 3
  t = 10
  t_list = get_t_list( t, h_list )
  x = exact_x_problem_1( t_list )
  dx = exact_xd_problem_1( t_list )
  for i in range(0, n):
    print( 'x({:.10e})={:.10e}'.format(t_list[i], x[i]) )
  for i in range(0, n):
    print( 'dx({:.10e})={:.10e}'.format(t_list[i], dx[i]) )
  c = discretize_bdf ( order, h_list )
  c_diff = discretize_bdf_diff ( order, h_list )

  tp = 10
  poly = interp_poly ( order, t_list, tp )
  xp = interp ( poly, x )
  print( 'xp({:.10e})={:.10e}'.format(tp, xp) )
  print( 'exact_x({:.10e})={:.10e}'.format(tp, exact_x_problem_1(tp)) )

  # this Lagrange coefficient should be the same as BDF diff coefficient
  poly_diff = interp_diff_poly ( order, t_list, tp )
  xdp = interp ( poly_diff, x )
  print( 'xdp({:.10e})={:.10e}'.format(tp, xdp) )
  print( 'exact_xd({:.10e})={:.10e}'.format(tp, exact_xd_problem_1(tp)) )

# main

test_problem_1()


