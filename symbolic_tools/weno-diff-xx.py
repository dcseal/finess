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

        u_stencil = [u0, u1, u2] and index = 0

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

def compute_poly_fit( u_stencil, x ):
    """Compute an expression for the polynomial approximation (in the variable
    x) to the given stencil.
    
    This polynomial fits points (x,p(x)) through

        ( 1, u[0] ), ( 2, u[1]), ( 3, u[2]), \dots
"""

    from sympy.polys.polyfuncs import interpolate
    return sympy.poly( interpolate( u_stencil, x ) )

def compute_h1_norm( u_stencil, indices ):
    """Compute the H1-norm of a given stencil.  This routine computes the
    integral

        beta = \sum_l dx**(2l-1) \int_{xl}^{xr} d^l/dx^l p(x)\, dx.

    that measures the smoothness of a polynomial that fits points

        ( x_0, u[0] ), (x_1, u[1]), (x_2, u[2]), \dots

    The parameters xl,xr = indices[:] are used to define the inteverval for
    integration.
"""

    from sympy.polys.polyfuncs import interpolate
    from sympy.abc import xi

    # Fit a polynomial through the whole stencil.
    # This fits the points, (x=1, p[0]), (x=2, p[1]), ...
    # p = compute_poly_fit( u_stencil, xi )
    p = sympy.poly( interpolate( u_stencil, xi ), xi )
    p = p.diff( xi )
    p = p.diff( xi )

#   print( ' ' )
#   print( u_stencil )
#   print( p.subs( xi, 1 ), p.subs( xi, 2 ), p.subs( xi, 3 ) )
#   print('p   = ', p   )
#   print('dp  = ', dp  )
#   print('d2p = ', d2p )
#   print(' ')

#   tmp = dpsqd.integrate( (xi, indices[0], indices[1] ) )
#   tmp = tmp + d2psqd.integrate( (xi, indices[0], indices[1] ) )
#   return tmp
    tmp = 0
    for mp in range( len( u_stencil ) ):
        tmp = tmp + (p**2).integrate( xi ).eval( xi, indices[1] ) - (p**2).integrate( xi ).eval( xi, indices[0] )
        p = p.diff( xi )
    return tmp

# (uniform) grid spacing
dx = sympy.symbols("dx")

uim3 = sympy.symbols("uim3")
uim2 = sympy.symbols("uim2")
uim1 = sympy.symbols("uim1")
ui   = sympy.symbols("ui")
uip1 = sympy.symbols("uip1")
uip2 = sympy.symbols("uip2")
uip3 = sympy.symbols("uip3")

# Three Lagrange polynomials and their derivatives:
beta0 = compute_h1_norm( [uim2, uim1, ui], (Rational(5,2), Rational(7,2) ) ) 
beta1 = compute_h1_norm( [uim1, ui, uip1], (Rational(3,2), Rational(5,2) ) ) 
beta2 = compute_h1_norm( [ui, uip1, uip2], (Rational(1,2), Rational(3,2) ) ) 

print('beta0 = ', beta0 )
print('beta1 = ', beta1 )
print('beta2 = ', beta2 )

# Exact smoothness indicators:
beta = [None]*3
beta[0] = Rational(13,12)*(uim2-2*uim1+ui)**2 + Rational(1,4)*(uim2-4*uim1+3*ui)**2
beta[1] = Rational(13,12)*(uim1-2*ui+uip1)**2 + Rational(1,4)*(uim1-uip1)**2
beta[2] = Rational(13,12)*(ui-2*uip1+uip2)**2 + Rational(1,4)*(3*ui-4*uip1+uip2)**2

print('exact beta[0] = ', beta[0])
print('exact beta[1] = ', beta[1])
print('exact beta[2] = ', beta[2])
print( sympy.simplify( sympy.expand( beta[0] - beta0 ) ) )

# Now, work out 2nd-derivative using larger stencil
u_stencil = [ uim3, uim2, uim1, ui, uip1, uip2, uip3 ]

# Compute derivative using the whole stencil:
gamma = compute_derivs( u_stencil, 2, dx )

# Four sub-stencils (of length four)
u0 = [uim3, uim2, uim1, ui  ]
u1 = [uim2, uim1,   ui, uip1]
u2 = [uim1, ui,   uip1, uip2]
u3 = [ui,   uip1, uip2, uip3]
u  = [u0, u1, u2, u3 ]

betax0 = compute_h1_norm( u0, (Rational(7,2), Rational(9,2) ) )
betax1 = compute_h1_norm( u1, (Rational(5,2), Rational(7,2) ) )
betax2 = compute_h1_norm( u2, (Rational(3,2), Rational(5,2) ) )
betax3 = compute_h1_norm( u3, (Rational(1,2), Rational(3,2) ) )

print('betax0 = ', betax0 )
print('betax1 = ', betax1 )
print('betax2 = ', betax2 )
print('betax3 = ', betax3 )


