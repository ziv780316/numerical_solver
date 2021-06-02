#!/usr/bin/env python3.8
import math
import numpy as np

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
  A[0, 0] = 1
  A[0, 1] = 1
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

  c = np.linalg.solve(A, b)

  if debug:
    for i in range(1, n):
      print( '\u03b1{id} = {c}'.format(id=i-1, c=c[i-1,0]), );
    print( '\u03b20 = {c}'.format(c=c[n-1,0]), );

  return c

# interpolation polynomial
def interp_poly ( order, t_list, t ):
  '''
  Lagrange Interpolation
  '''
  debug = True
  n = order + 1
  l = np.ones( (n,n) )
  for i in range(0, n):
    for j in range(0, n):
      if j == i:
        continue 
      l[i] *= (t - t_list[j])/(t_list[i]-t_list[j])

  if debug:
    for i in range(1, n):
      print( 'L{id} = {c}'.format(id=i-1, c=l[i-1,0]), );

  return l

def interp_diff_poly ( order, t_list, t ):
  '''
  Lagrange Differentiation Interpolation
  '''
  debug = True
  n = order + 1
  ld = np.ones( (n,n) )
  for i in range(0, n):
    numer = 0
    for j in range(0, n):
      if j == i:
        continue 
      ld[i] *= 1/(t_list[i]-t_list[j])
      tmp = 1
      for k in range(0, n):
        if k == i or k == j:
          continue 
        tmp *= (t - t_list[k])
      numer += tmp
    ld[i] *= numer

  if debug:
    for i in range(1, n):
      print( 'Ld{id} = {c}'.format(id=i-1, c=ld[i-1,0]), );

  return ld

# interpolation by polynomial
def interp ( l, x ):
  '''
  Lagrange Interpolation
  '''
  n = l.ndim
  xp = 0
  for i in range(0, n):
      xp += l[i,0]*x[i,0]

  return xp
  

# test problem
def test_f1 ( t, dt, x, xd ):
  lamda = 1
  f = (x + dt*lamda*xd) 
  return f

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


# main

order = 2
n = order + 1
h_list = np.zeros( (n,1) )
t_list = np.zeros( (n,1) )
x = np.zeros( (n,1) )
t = 3
h_list[0] = 1
h_list[1] = 1
h_list[2] = 1
x[0]=4
x[1]=1
x[2]=0
t_list[0] = t - h_list[0]
for i in range(1, n):
  t_list[i] = t_list[i-1] - h_list[i]
print( 't_list=\n', t_list )

discretize_bdf ( order, h_list )

tp = 4
poly = interp_poly ( order, t_list, tp )
poly_diff = interp_diff_poly ( order, t_list, tp )
xp = interp ( poly, x )
print( 'xp=', xp )
xdp = interp ( poly_diff, x )
print( 'xdp=', xdp )


