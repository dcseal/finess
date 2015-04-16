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

gamma5  = compute_derivs( np.arange(5 ), 2, dx )
gamma7  = compute_derivs( np.arange(7 ), 3, dx )
gamma9  = compute_derivs( np.arange(9 ), 4, dx )
gamma11 = compute_derivs( np.arange(11), 5, dx )
#gamma13 = compute_derivs( np.arange(13), 6, dx )

#sub_sten_val = g0*interps2[0] + g1*interps2[1] + g2*interps2[2]
#linear_val   = sum( (dx**2)*u*g for (u,g) in zip( u_stencil, gamma[2,:] ) )

stencil_size = 5
print("const double deriv_matrix%d[%d][%d] = {" % (stencil_size,stencil_size-1,stencil_size) )
for i in range(1,stencil_size):
    print("{ ")
    for j in range(stencil_size-1):
        print( '%2.25e, ' % float( gamma5[i,j]*dx**i ) )
    print( '%2.25e ' % float( gamma5[i,stencil_size-1]*dx**i ) )
    print('},')
print("};")

stencil_size = 7
print("const double deriv_matrix%d[%d][%d] = {" % (stencil_size,stencil_size-1,stencil_size) )
for i in range(1,stencil_size):
    print("{ ")
    for j in range(stencil_size-1):
        print( '%2.25e, ' % float( gamma7[i,j]*dx**i ) )
    print( '%2.25e ' % float( gamma7[i,stencil_size-1]*dx**i ) )
    print('},')
print("};")

stencil_size = 9
print("const double deriv_matrix%d[%d][%d] = {" % (stencil_size,stencil_size-1,stencil_size) )
for i in range(1,stencil_size):
    print("{ ")
    for j in range(stencil_size-1):
        print( '%2.25e, ' % float( gamma9[i,j]*dx**i ) )
    print( '%2.25e ' % float( gamma9[i,stencil_size-1]*dx**i ) )
    print('},')
print("};")

stencil_size = 11
print("const double deriv_matrix%d[%d][%d] = {" % (stencil_size,stencil_size-1,stencil_size) )
for i in range(1,stencil_size):
    print("{ ")
    for j in range(stencil_size-1):
        print( '%2.25e, ' % float( gamma11[i,j]*dx**i ) )
    print( '%2.25e ' % float( gamma11[i,stencil_size-1]*dx**i ) )
    print('},')
print("};")


