#==============================================================================#
# This file is part of HYPERPYWS: Hyperbolic Python WENO Solver
#
#   *** This software is made available "as is" without any assurance that it
#   *** will work for your purposes.  The software may in fact have defects, so
#   *** use the software at your own risk.
#
# License: GPL, see COPYING for details
#
# Copyright (C) 2013 
#
#    David Seal,  seal@math.msu.edu,  Michigan State University
#    Yaman Guclu, guclu@math.msu.edu, Michigan State University
#
#===============================================================================

def burger_exact(x, t, q0, qp):
    """ Solve for the exact solution of Burger's equation.
    
    Parameters:
    ===========

        x  : np array of points where we want the solution.

        t  : final time for the solution

        q0 : callable function describing initial conditions.  
            That is, q0(x) = q(t=0, x ).

        qp : callable function describing the spatial derivative of the initial 
             conditions.
    Returns:
    ========

        q  : exact solution evaluated at x, t.

    """

    # functions required for newton iteration:
    def f(z):
        return q0(z)*t + z - x
    def fp(z):
        return qp(z)*t + 1.0

    # newton convergence parameters
    MAX_ITER = 100
    tol      = 1e-15

    print('Computing exact solution')

    # initial guess:
    xi = x

    # Newton iteration until we find desired tolerance of EVERY point.
    num_iter = 0
    while( max( abs( f(xi) ) ) > tol ):
        xi = xi - ( f(xi) / fp(xi) )
        num_iter += 1
        if( num_iter > MAX_ITER ):
            print('  Too many iterations = %d ' % num_iter )
            break

    print(' Found exact solution in %d iterations' % num_iter )

    # exact solution:
    return q0(x-t*q0(xi))

