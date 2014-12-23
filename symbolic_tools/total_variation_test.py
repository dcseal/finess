"""Test total variation.

The purpose of this module is to test WENO methods on a single step function
example.  The WENO method has been implemented symbolically, and we wish to
study the total variation on a variety of methods.

In particular, the original motivation for studying these schemes was to look
at multiderivative formulations of the Picard integral formulation of the PDE.
"""

import sympy
from sympy import Rational

def reconstruct_left( eps, p, u_stencil ):
    """ Reconstruct u_{i+1/2}.
    """
    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # Compute smoothness indicators (identical for left/right values):
    beta = [None]*3
    beta[0]=sympy.Rational(13,12)*(uim2-2*uim1+ui)**2+sympy.Rational(1,4)*(uim2-4*uim1+3*ui)**2
    beta[1]=sympy.Rational(13,12)*(uim1-2*ui+uip1)**2+sympy.Rational(1,4)*(uim1-uip1)**2
    beta[2]=sympy.Rational(13,12)*(ui-2*uip1+uip2)**2+sympy.Rational(1,4)*(3*ui-4*uip1+uip2)**2
    
    # 3rd-order reconstructions using small 3-point stencils
    u1 = Rational( 1,3)*uim2 - Rational(7,6)*uim1 + Rational(11,6)*ui
    u2 = Rational(-1,6)*uim1 + Rational(5,6)*ui   + Rational( 1,3)*uip1
    u3 = Rational( 1,3)*ui   + Rational(5,6)*uip1 - Rational( 1,6)*uip2
    
    # Get linear weights and regularization parameter
    gamma = [Rational(1,10), Rational(6,10), Rational(3,10)]
    
    # Compute nonlinear weights and normalize their sum to 1
    omt  = [ g/(eps+b)**2 for g,b in zip(gamma,beta) ]
    omts = sum(omt)
    om   = [ o / omts for o in omt ]
    
    # Return 5th-order conservative reconstruction
    return om[0]*u1 + om[1]*u2 + om[2]*u3

def diff1NC( eps, p, u_stencil ):
    """ Compute a non-conservative difference using WENO weights.
    
    We compute u_x( x_i ) using a five-point stencil,

        {u_{i-2}, u_{i-1}, u_i, u_{i+1}, u_{i+2}

    This is a 'non-conservative' operator, because we do not attempt to
    reconstruct a primitive function from cell averages.
    """

    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # Compute smoothness indicators (identical for left/right values):
    beta = [None]*3
    beta[0]=sympy.Rational(13,12)*(uim2-2*uim1+ui)**2+sympy.Rational(1,4)*(uim2-4*uim1+3*ui)**2
    beta[1]=sympy.Rational(13,12)*(uim1-2*ui+uip1)**2+sympy.Rational(1,4)*(uim1-uip1)**2
    beta[2]=sympy.Rational(13,12)*(ui-2*uip1+uip2)**2+sympy.Rational(1,4)*(3*ui-4*uip1+uip2)**2

    # 3rd-order reconstructions using small 3-point stencils
    u1 = Rational( 1,2)*uim2 - 2*uim1 + Rational(3,2)*ui
    u2 = Rational(-1,2)*uim1          + Rational(1,2)*uip1
    u3 = Rational(-3,2)*ui   + 2*uip1 - Rational(1,2)*uip2
    
    # Get linear weights and regularization parameter
    gamma = [Rational(1,6), Rational(2,3), Rational(1,6)]
    
    # Compute nonlinear weights and normalize their sum to 1
    omt  = [ g/(eps+b)**2 for g,b in zip(gamma,beta) ]
    omts = sum(omt)
    om   = [ o / omts for o in omt ]
    
    # Return 5th-order conservative reconstruction
    return om[0]*u1 + om[1]*u2 + om[2]*u3

def diff1( eps, p, u_stencil ):
    """ Compute a centered difference of u on a five point stencil,
    
        {u_{i-2}, u_{i-1}, u_i, u_{i+1}, u_{i+2}

    We return an approximation to u_x( x_i ).
    """

    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # 3rd-order reconstructions using small 3-point stencils
    u1 = Rational( 1,2)*uim2 - 2*uim1 + Rational(3,2)*ui
    u2 = Rational(-1,2)*uim1          + Rational(1,2)*uip1
    u3 = Rational(-3,2)*ui   + 2*uip1 - Rational(1,2)*uip2
    
    # Get linear weights and regularization parameter
    gamma = [Rational(1,6), Rational(2,3), Rational(1,6)]
    
    # Return 5th-order conservative reconstruction
    return gamma[0]*u1 + gamma[1]*u2 + gamma[2]*u3


def single_step_small_stencil( eps, eps_num ):
    """Take a single Euler step on a small stencil with discontinuous initial
    data.  This function is deprecated in favor of ConstructL and EulerStep.
    """

    # Consider initial conditions, defined by
    #
    #      u_i = 1, i > 0, and u_i = 0, otherwise

    F0 = [0, 0, 0, 0, 0]        # i = -4, -3, -2, -1, 0, centered at i = -2
    F1 = [0, 0, 0, 0, 1]        # i = -3, -2, -1, 0,  1, centered at i = -1
    F2 = [0, 0, 0, 1, 1]        # i = -2, -1, 0, 1, 2,   centered at i =  0
    F3 = [0, 0, 1, 1, 1]        # i = -1, 0, 1, 2, 3,    centered at i =  1
    F4 = [0, 1, 1, 1, 1]        # i = 0, 1, 2, 3, 4,     centered at i =  2
    F5 = [1, 1, 1, 1, 1]        # i = 1, 2, 3, 4, 5,     centered at i =  3
    F6 = [1, 1, 1, 1, 1]        # i = 2, 3, 4, 5, 6,     centered at i =  4

    # Both of these are not needed, because we will be multiplying the result by
    # nu, the CFL number
    #dx  = sympy.symbols('dx')
    #dt  = sympy.symbols('dt')

    Fp0 = sympy.simplify( reconstruct_left(eps,2, F0 ) )
    Fp1 = sympy.simplify( reconstruct_left(eps,2, F1 ) )
    Fp2 = sympy.simplify( reconstruct_left(eps,2, F2 ) )
    Fp3 = sympy.simplify( reconstruct_left(eps,2, F3 ) )
    Fp4 = sympy.simplify( reconstruct_left(eps,2, F4 ) )
    Fp5 = sympy.simplify( reconstruct_left(eps,2, F5 ) )
    Fp6 = sympy.simplify( reconstruct_left(eps,2, F6 ) )
    Fp  = [Fp0, Fp1, Fp2, Fp3, Fp4 ]

    dx_times_Fx = sympy.simplify( [ (Fp1-Fp0), (Fp2-Fp1), (Fp3-Fp2), (Fp4-Fp3), (Fp5-Fp4), (Fp6-Fp5) ] )

    nu = sympy.symbols('nu')
    # i=-1,0,1,2,3,4
    Un   = [0, 0, 1, 1, 1, 1]

    Unp1Ap = []
    Unp1   = []
    for (i,u) in enumerate( Un ):
        Unp1Ap.append( sympy.simplify( u - nu * dx_times_Fx[ i ] ).subs(eps, eps_num))
        Unp1.append( sympy.simplify( u - nu * dx_times_Fx[ i ] ) )

    return Unp1, Unp1Ap

def ConstructL( Un, eps ):
    """Construct the right hand side of q_t = L(q) with the WENO method.

    We will do the entire thing symbolically.
    """

    import numpy as np

    mx    = len( Un )
    mbc   = 3

    # Flux function, with extra padding added in (zeroth-order extrapolation)
    F_Vec = np.concatenate( (np.array([0,0,0]), Un    ) )
    F_Vec = np.concatenate( (F_Vec, np.array([1,1,1]) ) )

    # Flux function.  Fp[0] = left-most flux value.  Fp[mx] = right-most value
    Fp = np.zeros( mx+1, dtype=object )
    for i in range(mx+1):
        
        # Pull the current stencil
        u_stencil = F_Vec[i:i+2*mbc-1]
        Fp[i] = sympy.simplify( reconstruct_left( eps, 2, u_stencil ) )

    # Flux values
    Fp = np.array( Fp )

    # Construct right-hand side values:
    L = np.zeros( mx, dtype=object )
    for i in range(mx):
        L[i] = sympy.simplify( -(Fp[i+1]-Fp[i]) )
    return L, Fp

def ConstructPIFL( Un, Unx, nu, eps ):
    """Construct right hand side for Picard Integral Formulation.

    This routine computs L = WENO( Un - dt*Un_x ).
    """

    import numpy as np

    mx    = len( Un )
    mbc   = 3

    # Flux function, with extra padding added in (zeroth-order extrapolation)
    F_Vec = np.concatenate( (np.array([0,0,0]), Un-nu*Unx/2 ) )
    F_Vec = np.concatenate( (F_Vec, np.array([1,1,1])       ) )

    # Flux function.  Fp[0] = left-most flux value.  Fp[mx] = right-most value
    Fp = np.zeros( mx+1, dtype=object )
    for i in range(mx+1):
        
        # Pull the current stencil
        u_stencil = F_Vec[i:i+2*mbc-1]
        Fp[i] = sympy.simplify( reconstruct_left( eps, 2, u_stencil ) )

    # Flux values
    Fp = np.array( Fp )

    # Construct right-hand side values:
    L = np.zeros( mx, dtype=object )
    for i in range(mx):
        L[i] = sympy.simplify( -(Fp[i+1]-Fp[i]) )
    return L, Fp

def ConstructUx( Un, eps ):
    """Compute a derivative of Un.
    """

    import numpy as np

    mx    = len( Un )
    mbc   = 3

    # Flux function, with extra padding added in (zeroth-order extrapolation)
    BigU = np.concatenate( (np.array([Un[0],Un[0]]), Un     ) )
    BigU = np.concatenate( (BigU, np.array([Un[-1],Un[-1]]) ) )

    Ux = np.zeros( mx, dtype=object )
    for i in range(mx):
        
        # Pull the current stencil
        u_stencil = BigU[i:i+2*mbc-1]
        Ux[i]     = sympy.simplify( diff1( eps, 2, u_stencil ) )

    return Ux 

def ConstructUxNC( Un, eps ):
    """Compute a derivative of Un.
    """

    import numpy as np

    mx    = len( Un )
    mbc   = 3

    # Flux function, with extra padding added in (zeroth-order extrapolation)
    BigU = np.concatenate( (np.array([Un[0],Un[0]]), Un     ) )
    BigU = np.concatenate( (BigU, np.array([Un[-1],Un[-1]]) ) )

    Ux = np.zeros( mx, dtype=object )
    for i in range(mx):
        
        # Pull the current stencil
        u_stencil = BigU[i:i+2*mbc-1]
        Ux[i]     = sympy.limit( sympy.simplify( diff1NC( eps, 2, u_stencil )), eps, 0 )
#       Ux[i]     = sympy.simplify( diff1NC( eps, 2, u_stencil ) )

    return Ux 


def EulerStep( Un, L, nu ):
    """Take a single Euler step."""

    Unp1   = np.zeros( mx, dtype=object )
    for i in range(mx):
        Unp1[i] = sympy.simplify( Un[i] + nu*L[i] )
    return Unp1

def TaylorStep( Un, Ut, Utt, nu ):
    """Take a single Euler step."""

    Unp1   = np.zeros( mx, dtype=object )
    for i in range(mx):
        Unp1[i] = sympy.simplify( Un[i] + nu*Ut[i] + nu**2/2*Utt[i] )
    return Unp1

nu      = sympy.symbols('nu')
#nu      = 1.0
eps     = sympy.symbols('eps')
#eps     = 1e-29

#Unp1Sm, Unp1ApSm = single_step_small_stencil( eps, eps_num )

###############################
# More general method
###############################

import numpy as np

# Consider initial conditions, defined by
#
#      u_i = 1, i > 0, and u_i = 0, otherwise
#
# Place half of the points to the left, and half to the right
mx   = 30
Un   = np.concatenate( (np.zeros( mx/2, dtype=int),   np.ones( mx/2, dtype=int   )) )
UnAp = np.concatenate( (np.zeros( mx/2, dtype=float), np.ones( mx/2, dtype=float )) )

# Construct right hand side for an Euler step

print("Constructing an Euler RHS")
L,Fp = ConstructL( Un, eps )
print("Taking an Euler step")
UFE = EulerStep( Un, L, nu )

Unx   = ConstructUx  ( Un, eps )
UnxNC = ConstructUxNC( Un, eps )

"""
print("Constructing a time-averaged RHS (CFD)")
Lt,Fpt      = ConstructPIFL( Un, Unx, nu, eps )
print("Taking an Euler step on time-averaged RHS (CFD)")
U_TaylorCFD = EulerStep( Un, Lt, nu )
"""

print("Constructing a time-averaged RHS (NCW)")
LtNC,FptNC  = ConstructPIFL( Un, UnxNC, nu, eps )
print("Taking an Euler step on time-averaged RHS (NCW)")
U_TaylorNCW = EulerStep( Un, LtNC, nu )

"""
print("Constructing Utt directly")
Utt = ConstructUx( -L, eps )
print("Taking Taylor step from Ut and Utt")
U_TaylorQS = TaylorStep( Un, L, Utt, nu )
"""

print("printing u from Forward-Euler time stepping")
for u in UFE:
    print( 'u = ', u, sympy.limit( u, eps, 0 ), ';' )

"""
print("printing U (CFD)")
for u in U_TaylorCFD:
    print( 'u = ', u, sympy.limit( u, eps, 0 ), ';' )
"""

print("printing U (NCW)")
for u in U_TaylorNCW:
    print( 'u = ', u, sympy.limit( u, eps, 0 ), ';' )

"""
print("printing U (QiuShu)")
for u in U_TaylorQS:
    print( 'u = ', u, sympy.limit( u, eps, 0 ), ';' )
"""

#print("Constructing new right hand side value")
#Lstar, Fpstar = ConstructL( Ustar, eps )

#Unp1 = Ustar
#Unp1  = Un + nu/2*( L + Lstar )

# Quick sanity check, do we get a reasonable value here!?
#print("printing Unp1")
#for u in Unp1:
#    print( u )
#   print( sympy.limit( u, eps, 0 ) )


#print("printing Ux")
#for i in range( len(UxNC) ):
#   print( ux )
#    print('UxNC = ', UxNC[i], 'WENO[u] = ', L[i] )
