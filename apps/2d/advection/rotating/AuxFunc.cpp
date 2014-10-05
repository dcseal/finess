#include "dogdefs.h"

// This is a user-required routine that defines the auxilary arrays when the
// parameter maux > 0.
//
// Each application is REQUIRED to define one of these.
//
//      Simple Advection equation.
//
// Input:
//
//    xpts( 1:numpts, 1:2 )   - The x, and y-coordinates for a list of points
//
// Output:
//
//    auxvals( 1:numpts, 1:maux )  - The vector containg auxilary values.
//
// See also: QinitFunc.
void AuxFunc(const dTensor2& xpts, dTensor2& auxvals)
{

    const int numpts = xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        // u:  1-component of the advection velocity
        auxvals.set(i,1,  2.0*pi*(y-0.5) );

        // v:  2-component of the advection velocity
        auxvals.set(i,2, -2.0*pi*(x-0.5) );

        // divu:  u_x + v_y
        auxvals.set(i,3,   0.0 );
    }

}
