#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
//      2D Euler equations.
//
//      Q(1) = density
//      Q(2) = rho u^1
//      Q(3) = rho u^2
//      Q(4) = rho u^3
//      Q(5) = Energy
//
// Input:
//
//    xpts( 1:numpts, 1:2 )      - The x, and y-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: $FINESS/lib/2d/blanks/QinitFunc
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

    // Gas consant
    const double gamma = global_ini_params.get_gamma();
    const double x0 = global_ini_params.get_x0();

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double xstar = x0 + y*osq3;

        double rho,u1,u2,u3,press;

        if(x < xstar)
        {
            rho   =  8.0;
            u1    =  8.25*cos(pi/6.0);
            u2    = -8.25*sin(pi/6.0);
            u3    =  0.0;
            press =  116.5;
        }
        else
        {
            rho   =  1.4;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            press =  1.0;
        } 

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
