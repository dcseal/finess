from __future__ import print_function  # For printing no newline
import sympy  # symbolic library
import re     # regular expressions: for saving useful output

def fixPowers(s):
    return re.sub(r'q(\d+)\*\*(\d+)', r'pow( q\1, \2 )', s)

meqn = 3

# Conserved variables:
q1 = sympy.symbols("q1")
q2 = sympy.symbols("q2")
q3 = sympy.symbols("q3")

# Flux function -- first component
f1 = sympy.symbols("f1")
f2 = sympy.symbols("f2")
f3 = sympy.symbols("f3")

# Flux function -- second component
g1 = sympy.symbols("g1")
g2 = sympy.symbols("g2")
g3 = sympy.symbols("g3")


# Primitive variables:
h  = sympy.symbols("h")
u1 = sympy.symbols("u1")
u2 = sympy.symbols("u2")

h  = q1
u1 = q2/h
u2 = q3/h

f1 = h*u1
f2 = h*u1**2 + (h**2)/2
f3 = h*u1*u2

g1 = h*u2
g2 = h*u1*u2
g3 = h*u2**2 + (h**2)/2

Q = [q1, q2, q3]
F = [f1, f2, f3]
G = [g1, g2, g3]

# A = sympy.Matrix( meqn, meqn )
print("Computing the Jacobian of the flux function, f'(q)")
for j in range(meqn):
  for k in range(meqn):
      print( ('Dflux.set(i,%d,%d,1, ' % (j+1,k+1) ), end="" )
      tmp = fixPowers( str(sympy.simplify( sympy.expand( sympy.diff( F[j], Q[k]))) ) )
      print( tmp, end=");\n")
print(' ')

print("Computing the Jacobian of the flux function, g'(q)")
for j in range(meqn):
  for k in range(meqn):
      print( ('Dflux.set(i,%d,%d,2, ' % (j+1,k+1) ), end="" )
      tmp = fixPowers( str(sympy.simplify( sympy.expand( sympy.diff( G[j], Q[k]))) ) )
      print( tmp, end=");\n")
print(' ')

print("Computing the Hessian of the flux function: f''(q)")
for m1 in range(meqn):
  print(' ')
  for m2 in range(meqn):
    for m3 in range(meqn):
      print( ('D2flux.set(i,%d,%d,%d,1, ' % (m1+1,m2+1,m3+1) ), end="" )
      tmp = fixPowers( str( sympy.expand( sympy.diff( F[m1], Q[m2], Q[m3]))) )
      print( tmp, end=");\n")
print(' ')

print("Computing the Hessian of the flux function: g''(q)")
for m1 in range(meqn):
  print(' ')
  for m2 in range(meqn):
    for m3 in range(meqn):
      print( ('D2flux.set(i,%d,%d,%d,2, ' % (m1+1,m2+1,m3+1) ), end="" )
      tmp = fixPowers( str( sympy.expand( sympy.diff( G[m1], Q[m2], Q[m3]))) )
      print( tmp, end=");\n")
print(' ')

