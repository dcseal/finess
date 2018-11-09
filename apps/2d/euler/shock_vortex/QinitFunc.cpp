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


// See 5.7. Shock–vortex interaction in "A posteriori subcell limiting of the
// discontinuous Galerkin finite element method for hyperbolic conservation
// laws," Dumbser et al, JCP 2014.

    // Gas consant
    const double gamma = global_ini_params.get_gamma();

    // Parameters for this probelm
    const double a    = 0.075; 
    const double b    = 0.175; 
    const double MS   = 1.5; 
    const double MV   = 0.7; 
    const double p0   = 1.; 
    const double rho0 = 1.;

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        const double xc = 0.25;
        const double yc = 0.5;

        double rho,u1,u2,u3,press;

        // Distance to center
        const double r = sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) );

        double vphi;

        // Check distance to center of vortex
        if( r < a )
        {
            vphi = vm * r / a;
        }
        else if( r < b )
        {
            vphi = vm * a / (a*a - b*b) * (r - b*b/r);
        } 
        else
        {
            vphi = 0.;
        }

        // Comptue derived quantities and store each of the variables
        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
