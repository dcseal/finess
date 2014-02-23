from __future__ import print_function  # For printing no newline
import sympy
import re
 
def fixPowers(s):
    return re.sub(r'q(\d+)\*\*(\d+)', r'pow( q\1, \2 )', s)

meqn = 5

# Ratio of specific heats
gamma = sympy.symbols("gamma")

q1 = sympy.symbols("q1")
q2 = sympy.symbols("q2")
q3 = sympy.symbols("q3")
q4 = sympy.symbols("q4")
q5 = sympy.symbols("q5")

f1 = sympy.symbols("f1")
f2 = sympy.symbols("f2")
f3 = sympy.symbols("f3")
f4 = sympy.symbols("f4")
f5 = sympy.symbols("f5")

g1 = sympy.symbols("g1")
g2 = sympy.symbols("g2")
g3 = sympy.symbols("g3")
g4 = sympy.symbols("g4")
g5 = sympy.symbols("g5")


# Primitive variables:
u1 = sympy.symbols("u1")
u2 = sympy.symbols("u2")
u3 = sympy.symbols("u3")
en = sympy.symbols("en")  # energy
pr = sympy.symbols("pr")  # pressure

u1 = q2/q1
u2 = q3/q1
u3 = q4/q1
pr = (gamma-1)*(en-q1/2*(u1**2+u2**2+u3**2))
en = q5

f1 = q1*u1
f2 = q1*u1**2 + pr
f3 = q1*u1*u2
f4 = q1*u1*u3
f5 = u1*(en+pr)

g1 = q1*u2
g2 = q1*u2*u1
g3 = q1*u2**2 + pr
g4 = q1*u2*u3
g5 = u2*(en+pr)

Q = [q1, q2, q3, q4, q5]
F = [f1, f2, f3, f4, f5]
G = [g1, g2, g3, g4, g5]

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
