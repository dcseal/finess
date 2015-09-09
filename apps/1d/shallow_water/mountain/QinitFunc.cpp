#include <cmath>
#include "tensors.h"
#include "stdio.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize();

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i);

printf("x = %f\n", x );
        // Bottom topography
        double bot = 0.5e0*exp(-100.0e0*pow((x-0.5),2));

        //h = 1.0040014495278291 - bot;
        //Qflow = 0.2;

        double h     = 1.0678715170821127 - bot;
        double Qflow = 0.25;

        qvals.set(i,1, h );      // depth
        qvals.set(i,2, Qflow );  // momentum

    }
}
