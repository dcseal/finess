#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// For Sedov blast problem, this is supposed to be a delta function in energy
// at x=y=0.  To approximate this, we simply set a single spike at the lower
// left most corner of the domain.
//
// This additional energy will be inserted in the routine, AfterQinit.
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
//  const double gamma = global_ini_params.get_gamma();
//  const double dx    = global_ini_params.get_dx();
//  const double dy    = global_ini_params.get_dy();

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        const double rho = 1.0;
        const double u1  = 0.0;
        const double u2  = 0.0;
        const double u3  = 0.0;
        const double energy = 1.0e-12;

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
