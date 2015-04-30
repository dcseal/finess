from __future__ import print_function  # For printing no newline
import sympy
from sympy import Rational
from sympy import factorial
import numpy as np
 
def Taylor( n, dx ):
    """Compyute n terms in the Taylor expansion for a function centered 'dx'
    away from where the terms will be evaluated."""

    return [ (dx)**j/factorial(j) for j in range(n) ]


def Taylor2d( n, I,J,dx,dy ):
    """Compyute n terms in the Taylor expansion for a function centered 'dx'
    away from where the terms will be evaluated."""
    v=np.zeros(I.flatten().shape)
    I1=I.flatten()
    J1=J.flatten()
    nx=I.shape[0]
    ny=I.shape[1]
    N=nx*ny
    #print("here1",nx)
    #print("here2",(dx)**I1[nx-1]*(dy)**J1[nx-1])    
    return [ (dx)**I1[j]/factorial(I1[j])*(dy)**J1[j]/factorial(J1[j]) for j in range(N) ]



def compute_derivs( n, index, dx, dy ):
    """Compute all of the derivatives of u_stencil.  Index refers to the
    location where these values will be computed.  For example,

        u_stencil = [u0, u1, u0] and index = 0

    will compute all of the derivatives at u0, u0_x, u0_xx.

    However, index = 0.5 will compute the derivatives u1/2, u1/2_x and
    u1/2_xx.
    """

    # Set up matrix that needs to be inverted
    N1=n*n;
    i1=np.array(range(n))
    [I,J]=np.meshgrid(i1,i1)
    I1=I.flatten()
    J1=J.flatten()

    A1 = sympy.Matrix( N1, N1, lambda i,j : 0 )
    A2 = sympy.Matrix( N1, N1, lambda i,j : 0 )

    Z = sympy.Matrix( N1, 2, lambda i,j : 0 )
    
    B=sympy.Matrix(N1,1, lambda i,j : 0);

    for i in range(N1):
            
            t = Taylor2d( n, I,J ,(I1[i]-index)*dx,(J1[i]-index)*dy)
            t1 = Taylor2d( n, I,J ,(I1[i]-index),(J1[i]-index))
            Z[i,0]=(I1[i]-index)*dx
            Z[i,1]=(J1[i]-index)*dy
            B[i]=(I1[i]-index)**3*(J1[i]-index)**2
            for k in range(N1):
                A1[i,k] = t[k]
                A2[i,k] = t1[k]
  
    A2=A2.inv()
    return A2

np.set_printoptions(precision=25)

# (uniform) grid spacing
dx = sympy.symbols("dx")
dy = sympy.symbols("dy")

A2 = compute_derivs( 7,3, dx, dy )
size1=A2.shape[0]
f1=open('diffMatrix7.h', 'w+')
print("const double Mdderiv_matrix7[%d][%d] = {" % (size1,size1),file=f1 )
for i in range(0,size1):
       print("{ ",file=f1)
       for j in range(size1-1):
        print( '%2.25e, ' % float( A2[i,j] ),file=f1 )
       print( '%2.25e ' % float( A2[i,size1-1] ),file=f1 )
       print('},',file=f1)
print("};",file=f1)

