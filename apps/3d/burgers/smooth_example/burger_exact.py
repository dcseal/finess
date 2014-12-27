# -*- coding: utf-8 -*-


def burger_exact_3d(x, y, z, t, q0, max_abs_sum_partial_deriv_q0):
    """ Solve for the exact solution of Burger's equation.
    
    Parameters:
    ===========

        x  : np array of points where we want the solution.

        t  : final time for the solution

        q0 : callable function describing initial conditions.  
            That is, q0(x, y, z) = q(t=0, x, y, z ).

        max_abs_sum_partial_deriv_q0:
            max of abs(\partial q0 / \partial x + ...)
            Used for estimating convergence.  Can overestimate a
            little bit.

    Returns:
    ========

        q  : exact solution evaluated at x, y, z, t.

    """
    assert t >= 0
    # Tolerance for error.
    tol      = 1e-15

    from math import ceil, log

    a = max_abs_sum_partial_deriv_q0 * t
    assert a < 1

    # function required for fixed point iteration:
    def f(q):
        return q0(x - q*t, y - q*t, z - q*t)
    
    qzeroth = q0(x, y, z)
    qfirst  = f(qzeroth)
    d1 = abs(qfirst - qzeroth)
    if d1 < tol:
        return qzeroth

    NUM_ITER = int(ceil(log((1.0-a) * tol / d1) / log(a)))

    q = qzeroth
    for i in range(NUM_ITER):
        q = f(q)

    return q
#
#
#from pylab import linspace, plot, empty_like
#from math import sin, pi
#
#X = linspace(0, 2, num = 5000)
#q0 = lambda x: 0.5 + sin(pi * x)
#max_magnitude_of_divergence_of_q0 = pi
#t = 1 / pi * 0.5
#
#exact_solution = empty_like(X)
#
#for i in range(X.size):
#    exact_solution[i] = burger_exact_1d(X[i], t, q0, max_magnitude_of_divergence_of_q0)
#
#plot(X, exact_solution)
