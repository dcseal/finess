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

eps     = sympy.symbols('eps')
eps_num = 1e-29

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
Unp1   = []
Unp1Ap = []
for (i,u) in enumerate( Un ):
    Unp1Ap.append( sympy.simplify( u - nu * dx_times_Fx[ i ] ).subs(eps, eps_num))
    Unp1.append( sympy.simplify( u - nu * dx_times_Fx[ i ] ) )

#def compute_tv_euler_step( nu ):
tv = 0
#for i in range(1,len(Unp1)):
#    tv = tv + abs( Unp1[i].subs( eps, 1e-29 ) - Unp1[i-1].subs( eps, 1e-29 ) )

#print( tv )




