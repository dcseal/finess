"""
The purpose of this module is to provide access to the usual WENO and difference 
operators found throghout FINESS.

Required Modules:
-----------------

    * sympy - symbolic toolbox for Python.
"""

import sympy

def central_diff1( u, dx ):
    """Compute first-deriative using a five-point centered stencil."""
    u_x  = (  u[0] - 8*(u[1]-u[3]) - u[4] )/(12*dx)
    return u_x

def central_diff2( u, dx ):
    """Compute second-deriative using a five-point centered stencil."""
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

def LinearReconstruct( u_stencil ):
    """WENO reconstruction with linear weights.  This reconstructs u_{i+1/2}
    given cell averages \\bar{u}_i at each neighboring location.
    
    See `High Order Weighted Essentially Nonoscillatory Schemes for
    Convection Dominated Problems'. 

    Input
    -----

        u_stencil : 
            stencil describing current solution

    Returns
    -------

        uiphf     :
            u_{i+1/2} after performing the full reconstruction procedure.

    See also: WENOReconstruct.
    """

    uim2, uim1, ui, uip1, uip2 = u_stencil
    return uim2/30 - (13*uim1)/60 + 47*(ui/60) + 9*(uip1/20) - uip2/20

def WENOReconstruct( u_stencil, eps=1.0e-29, p=2 ):
    """WENO reconstruction.  This reconstructs u_{i+1/2}
    given cell averages \\bar{u}_i at each neighboring location.
    
    See `High Order Weighted Essentially Nonoscillatory Schemes for
    Convection Dominated Problems'. 

    Input
    -----

        u_stencil : 
            stencil describing current solution
        eps       : 
            regularization parameter
        p         :
            power parameter

    Returns
    -------

        uiphf     :
            u_{i+1/2} after performing the full reconstruction procedure.

    """

    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # Compute smoothness indicators (identical for left/right values):
    beta = [None]*3
    beta[0]=(13./12.)*(uim2-2*uim1+ui)**2+0.25*(uim2-4*uim1+3*ui)**2
    beta[1]=(13./12.)*(uim1-2*ui+uip1)**2+0.25*(uim1-uip1)**2
    beta[2]=(13./12.)*(ui-2*uip1+uip2)**2+0.25*(3*ui-4*uip1+uip2)**2
    
    # 3rd-order reconstructions using small 3-point stencils
    u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui
    u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1
    u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2
    
    # Get linear weights and regularization parameter
    gamma = [0.1, 0.6, 0.3]
    
    # Compute nonlinear weights and normalize their sum to 1
    omt  = [ g/(eps+b)**2 for g,b in zip(gamma,beta) ]
    omts = sum(omt)
    om   = [ o / omts for o in omt ]
    
    # Return 5th-order conservative reconstruction
    return om[0]*u1 + om[1]*u2 + om[2]*u3
  
def WENODiffX( u_stencil, eps=1.0e-29, p=2 ):
    """Non-conservative WENO differentiation.  This defines u'_{i+1/2}
    given cell averages \\bar{u}_i at each neighboring location.
    
    See `High Order Weighted Essentially Nonoscillatory Schemes for
    Convection Dominated Problems'. 

    Input
    -----

        u_stencil : 
            stencil describing current solution
        eps       : 
            regularization parameter
        p         :
            power parameter

    Returns
    -------

        uiphf     :
            u_{i+1/2} after performing the full reconstruction procedure.

    """

    # lazy, but readable indexing into list that's being passed in:
    uim2, uim1, ui, uip1, uip2 = u_stencil
    
    # Compute smoothness indicators (identical for left/right values):
    beta = [None]*3
    beta[0]=(13./12.)*(uim2-2*uim1+ui)**2+0.25*(uim2-4*uim1+3*ui)**2
    beta[1]=(13./12.)*(uim1-2*ui+uip1)**2+0.25*(uim1-uip1)**2
    beta[2]=(13./12.)*(ui-2*uip1+uip2)**2+0.25*(3*ui-4*uip1+uip2)**2
    
    # 3rd-order reconstructions using small 3-point stencils
    u1 = ( 0.5  )*uim2 - (2.   )*uim1 + ( 1.5  )*ui;
    u2 = (-0.5  )*uim1 + (0.   )*ui   + ( 0.5  )*uip1;
    u3 = (-1.5  )*ui   + (2.   )*uip1 - ( 0.5  )*uip2;
    
    # Get linear weights and regularization parameter
    gamma = [1./6., 2./3., 1./6.]
    
    # Compute nonlinear weights and normalize their sum to 1
    omt  = [ g/(eps+b)**2 for g,b in zip(gamma,beta) ]
    omts = sum(omt)
    om   = [ o / omts for o in omt ]
    
    # Return 5th-order conservative reconstruction
    return om[0]*u1 + om[1]*u2 + om[2]*u3
  

