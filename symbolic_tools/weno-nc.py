from __future__ import print_function  # For printing no newline
import sympy
from sympy import Rational
from sympy import factorial
import numpy as np
 
def Taylor( n, dx ):
    """Compyute n terms in the Taylor expansion for a function centered 'dx'
    away from where the terms will be evaluated."""

    return [ (dx)**j/factorial(j) for j in range(n) ]

def compute_derivs( u_stencil, index, dx ):
    """Compute all of the derivatives of u_stencil.  Index refers to the
    location where these values will be computed.  For example,

        u_stencil = [u0, u1, u0] and index = 0

    will compute all of the derivatives at u0, u0_x, u0_xx.

    However, index = 0.5 will compute the derivatives u1/2, u1/2_x and
    u1/2_xx.
"""

    # Set up matrix that needs to be inverted
    n = len( u_stencil )
    A = sympy.Matrix( n, n, lambda i,j : 0 )
    for i in range(n):
        t = Taylor( n, (i-index)*dx )
        for j in range(n):
            A[i,j] = t[j]

    return A.inv()

# (uniform) grid spacing
dx = sympy.symbols("dx")

uim2 = sympy.symbols("uim2")
uim1 = sympy.symbols("uim1")
ui   = sympy.symbols("ui")
uip1 = sympy.symbols("uip1")
uip2 = sympy.symbols("uip2")

u_stencil = [uim2, uim1, ui, uip1, uip2]

# Compute derivative using the whole stencil:
gamma = compute_derivs( u_stencil, 2, dx )

# Three sub-stencils
u0 = [uim2, uim1,   ui]
u1 = [uim1, ui,   uip1]
u2 = [ui,   uip1, uip2]
u  = [u0, u1, u2]

# Three Lagrange polynomials and their derivatives:
p0 = compute_derivs( u0, 2, dx )
p1 = compute_derivs( u1, 1, dx )
p2 = compute_derivs( u2, 0, dx )
p  = [p0, p1, p2]

interps = [0,0,0]
for j in range(3):
    interps[j] = u[j][0]*(p[j][1,0]*dx) +  u[j][1]*(p[j][1,1]*dx) + u[j][2]*(p[j][1,2]*dx)

# Goal: find unknowns (g0, g1, g2) such that
#       g0 * u0 + g1 * u1 + g2 * u2 = gamma[1,:]*u_stencil
g0 = Rational(1,6)
g1 = Rational(2,3)
g2 = Rational(1,6)

sub_sten_val = g0*interps[0] + g1*interps[1] + g2*interps[2]
linear_val   = sum( dx*u*g for (u,g) in zip( u_stencil, gamma[1,:] ) )

error = sympy.simplify( sympy.expand( sub_sten_val - linear_val ) )
print( error )

interps2 = [0,0,0]
for j in range(3):
    interps2[j] = u[j][0]*(p[j][2,0]*dx**2) +  u[j][1]*(p[j][2,1]*dx**2) + u[j][2]*(p[j][2,2]*dx**2)

# Goal: find unknowns (g0, g1, g2) such that
#       g0 * u0 + g1 * u1 + g2 * u2 = gamma[1,:]*u_stencil
g0 = Rational(-1,12)
g1 = Rational( 7, 6)
g2 = Rational(-1,12)

sub_sten_val = g0*interps2[0] + g1*interps2[1] + g2*interps2[2]
linear_val   = sum( (dx**2)*u*g for (u,g) in zip( u_stencil, gamma[2,:] ) )

error = sympy.simplify( sympy.expand( sub_sten_val - linear_val ) )
print( error )

# Smoothness indicators
#   eps  = sympy.symbols("epsilon")
#   beta = [None]*3
#   beta[0] = Rational(13,12)*(uim2-2*uim1+ui)**2 + Rational(1,4)*(uim2-4*uim1+3*ui)**2
#   beta[1] = Rational(13,12)*(uim1-2*ui+uip1)**2 + Rational(1,4)*(uim1-uip1)**2
#   beta[2] = Rational(13,12)*(ui-2*uip1+uip2)**2 + Rational(1,4)*(3*ui-4*uip1+uip2)**2
