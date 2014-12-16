"""
The purpose of this module is to perform a single Euler step on a step
function with the WENO method.

The initial conditions are

    q(t=0, x) = 1_{ x > 0 }

That is, q = 1 if x > 0, and q = 0 otherwise.

Required Modules:
-----------------

    * sympy - symbolic toolbox for Python.
"""

from weno_diff import central_diff1, central_diff2, ConstructIntegratedF, LinearReconstruct

import sympy

# Quick accessors

dx = sympy.Symbol("dx", real=True )    # spatial resolution
dt = sympy.Symbol("dt", real=True )    # time step
nu = sympy.Symbol("nu", real=True )    # CFL number
u  = sympy.Symbol("u" , real=True )    # advection speed

dt = nu*dx/u

qim5 = 0
qim4 = 0
qim3 = 0
qim2 = 0
qim1 = 0
qi   = 1
qip1 = 1
qip2 = 1
qip3 = 1
qip4 = 1

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
#forward_euler = sympy.collect( 1 - (u*dt)*( qi - qim1 ) / dx, exp(-I*k*dx) )
#print('Stability polynomial for Forward Euler is ')
#sympy.pretty_print( forward_euler )

