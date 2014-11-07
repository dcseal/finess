from __future__ import print_function  # For printing no newline
import sympy
import re
 
def fixPowers(s):
    return re.sub(r'q(\d+)\*\*(\d+)', r'pow( q\1, \2 )', s)

meqn = 8


# Ratio of specific heats
gamma = sympy.symbols("gamma")

q1 = sympy.symbols("q1")
q2 = sympy.symbols("q2")
q3 = sympy.symbols("q3")
q4 = sympy.symbols("q4")
q5 = sympy.symbols("q5")
q6 = sympy.symbols("q6")
q7 = sympy.symbols("q7")
q8 = sympy.symbols("q8")

f1 = sympy.symbols("f1")
f2 = sympy.symbols("f2")
f3 = sympy.symbols("f3")
f4 = sympy.symbols("f4")
f5 = sympy.symbols("f5")
f6 = sympy.symbols("f6")
f7 = sympy.symbols("f7")
f8 = sympy.symbols("f8")

g1 = sympy.symbols("g1")
g2 = sympy.symbols("g2")
g3 = sympy.symbols("g3")
g4 = sympy.symbols("g4")
g5 = sympy.symbols("g5")
g6 = sympy.symbols("g6")
g7 = sympy.symbols("g7")
g8 = sympy.symbols("g8")


# Primitive variables:
u1 = sympy.symbols("u1")
u2 = sympy.symbols("u2")
u3 = sympy.symbols("u3")

rho = q1
u1 = q2/q1
u2 = q3/q1
u3 = q4/q1
energy = q5
B1 = q6
B2 = q7
B3 = q8

pressure = (gamma-1)*(energy - rho/2*(u1**2+u2**2+u3**2) - (B1*B1 + B2*B2 + B3*B3)/2)

Bm = (B1*B1 + B2*B2 + B3*B3) / 2
Bu = u1*B1 + u2*B2 + u3*B3

f1 = rho*u1
f2 = rho*u1**2 + pressure + Bm - B1*B1
f3 = rho*u1*u2 - B1*B2
f4 = rho*u1*u3 - B1*B3
f5 = u1*(energy + pressure + Bm) - B1*Bu
f6 = 0
f7 = u1*B2 - u2*B1
f8 = u1*B3 - u3*B1

g1 = rho*u2
g2 = rho*u1*u2 - B2*B1
g3 = rho*u2**2 + pressure + Bm - B2*B2
g4 = rho*u2*u3 - B2*B3
g5 = u2*(energy + pressure + Bm) - B2*Bu
g6 = u2*B1 - u1*B2
g7 = 0
g8 = u2*B3 - u3*B2

Q = [q1, q2, q3, q4, q5, q6, q7, q8]
F = [f1, f2, f3, f4, f5, f6, f7, f8]
G = [g1, g2, g3, g4, g5, g6, g7, g8]

# A = sympy.Matrix( meqn, meqn )
print("// Computing the Jacobian of the flux function, f'(q)")
for j in range(meqn):
  for k in range(meqn):
      print( ('Dflux.set(i,%d,%d,1, ' % (j+1,k+1) ), end="" )
      tmp = fixPowers( str(sympy.simplify( sympy.expand( sympy.diff( F[j], Q[k]))) ) )
      print( tmp, end=");\n")
print(' ')

print("// Computing the Jacobian of the flux function, g'(q)")
for j in range(meqn):
  for k in range(meqn):
      print( ('Dflux.set(i,%d,%d,2, ' % (j+1,k+1) ), end="" )
      tmp = fixPowers( str(sympy.simplify( sympy.expand( sympy.diff( G[j], Q[k]))) ) )
      print( tmp, end=");\n")
print(' ')

print("// Computing the Hessian of the flux function: f''(q)")
for m1 in range(meqn):
  print(' ')
  for m2 in range(meqn):
    for m3 in range(meqn):
      print( ('D2flux.set(i,%d,%d,%d,1, ' % (m1+1,m2+1,m3+1) ), end="" )
      tmp = fixPowers( str( sympy.expand( sympy.diff( F[m1], Q[m2], Q[m3]))) )
      print( tmp, end=");\n")
print(' ')

print("// Computing the Hessian of the flux function: g''(q)")
for m1 in range(meqn):
  print(' ')
  for m2 in range(meqn):
    for m3 in range(meqn):
      print( ('D2flux.set(i,%d,%d,%d,2, ' % (m1+1,m2+1,m3+1) ), end="" )
      tmp = fixPowers( str( sympy.expand( sympy.diff( G[m1], Q[m2], Q[m3]))) )
      print( tmp, end=");\n")
print(' ')
