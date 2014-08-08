"""
The purpose of this module is to perform a Von-Neumann stability analysis on
the proposed Lax-Wendroff formulation of the Picard Integral Formulation of
WENO.

Required Modules:
-----------------

    * sympy - symbolic toolbox for Python.
"""

import sympy

# Quick accessors
I   = sympy.I    # equals \sqrt{-1}
exp = sympy.exp  # exponential operator
pi  = sympy.pi   # mathematical pi

def central_diff1( u, dx ):
    """Compute first-deriative using a five-point centered stencil."""
    u_x  = (  u[0] - 8*(u[1]-u[3]) - u[4] )/(12*dx)
    return u_x

def central_diff2( u, dx ):
    """Compute second-deriative using a five-point centered stencil."""
#   u_xx = ( -u[4] + 16*( u[3] + u[1] ) - 30*u[2] - u[0] )/( 12*(dx**2) )
    u_xx = ( -u[4] + 16*( u[3] + u[1] ) - 30*u[2] - u[0] )/( 12*(dx**2) )
    return u_xx

def ConstructIntegratedF( q_stencil, u, dx, dt ):
    """Construct the integrated flux function at a single point

        F = f - dt/2 f_t + dt**2/6 f_tt.

    For a linear problem, f = uq, f_t = -u**2 q_x, and f_tt = u**3 q_{xx},
    which reduces to

        F = u*( q - (u*dt)/2 q_x + (u*dt)**2 / 6 * q_xx )

    We assume that a 5-point stencil is passed in, and this routine computes F
    at the single point located at the center of this stencil.
    """

    assert( len( q_stencil ) == 5 )
    q_x  = sympy.simplify( central_diff1( q_stencil, dx ) )
    q_xx = sympy.simplify( central_diff2( q_stencil, dx ) )

    F = u*( q_stencil[2] - (u*dt)*(q_x/2) + (u*dt)**2*(q_xx/6) )

    return sympy.simplify( F )

def LinearReconstruct( u_stencil ):
    """Linear reconstruction with linear weights.  This reconstructs u_{i+1/2}
    given cell averages \\bar{u}_i at each neighboring location.
    
    See (2.14) in `High Order Weighted Essentially Nonoscillatory Schemes for
    Convection Dominated Problems'. 
    """

    uim2, uim1, ui, uip1, uip2 = u_stencil

    return uim2/30 - (13*uim1)/60 + 47*(ui/60) + 9*(uip1/20) - uip2/20

dx = sympy.Symbol("dx", real=True )    # spatial resolution
dt = sympy.Symbol("dt", real=True )    # time step
nu = sympy.Symbol("nu", real=True )    # CFL number
u  = sympy.Symbol("u" , real=True )    # advection speed

dt = nu*dx/u

# Ansatz: u_i = exp( I * k * x ).  Note that x = j*dx for some index i.
# WLOG, we'll assume that the stencil is centered at i=0.
k  = sympy.symbols("k", real=True )  # Fourier number

qim5 = exp( I * k * (-5) * dx )
qim4 = exp( I * k * (-4) * dx )
qim3 = exp( I * k * (-3) * dx )
qim2 = exp( I * k * (-2) * dx )
qim1 = exp( I * k * (-1) * dx )
qi   = exp( I * k * ( 0) * dx )
qip1 = exp( I * k * ( 1) * dx )
qip2 = exp( I * k * ( 2) * dx )
qip3 = exp( I * k * ( 3) * dx )
qip4 = exp( I * k * ( 4) * dx )

# Expansion of f:
#
#    F \approx f_i - (u*dt)/2 u_{i,x} + (u*dt)^2/3 u_{i,xx}
#
Fim3 = sympy.simplify( ConstructIntegratedF( [ qim5, qim4, qim3, qim2, qim1 ], u, dx, dt ) )
Fim2 = sympy.simplify( ConstructIntegratedF( [ qim4, qim3, qim2, qim1, qi   ], u, dx, dt ) )
Fim1 = sympy.simplify( ConstructIntegratedF( [ qim3, qim2, qim1, qi  , qip1 ], u, dx, dt ) )
Fi   = sympy.simplify( ConstructIntegratedF( [ qim2, qim1, qi  , qip1, qip2 ], u, dx, dt ) )
Fip1 = sympy.simplify( ConstructIntegratedF( [ qim1, qi  , qip1, qip2, qip3 ], u, dx, dt ) )
Fip2 = sympy.simplify( ConstructIntegratedF( [ qi  , qip1, qip2, qip3, qip4 ], u, dx, dt ) )

Fimh = sympy.simplify( LinearReconstruct( [Fim3, Fim2, Fim1, Fi, Fip1] ) )
Fiph = sympy.simplify( LinearReconstruct( [Fim2, Fim1, Fi, Fip1, Fip2] ) )

# Update is q_i^{n+1} = q_i^n - dt ( Fiph - Fim2 ) / dx
stability_polynomial = sympy.collect( sympy.expand(1 - dt*(Fiph - Fimh)/dx), nu )
print('Stability polynomial for Taylor PIF-WENO is ')
sympy.pretty_print( stability_polynomial )
print('Stability polynomial for Taylor PIF-WENO is ')
print( stability_polynomial )

# Try this with a simpler scheme that has a known solution:
forward_euler = sympy.collect( 1 - (u*dt)*( qi - qim1 ) / dx, exp(-I*k*dx) )
print('Stability polynomial for Forward Euler is ')
sympy.pretty_print( forward_euler )

