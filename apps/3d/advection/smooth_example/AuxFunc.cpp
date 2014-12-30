#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor2& xpts, dTensor2& auxvals)
{

    const int numpts=xpts.getsize(1);

    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i,1);
        double y = xpts.get(i,2);
        double z = xpts.get(i,3);

        // u1:  1-component of the advection velocity
        auxvals.set(i,1, 1.0 );

        // u2:  2-component of the advection velocity
        auxvals.set(i,2, 1.0);

        // u3:  3-component of the advection velocity
        auxvals.set(i,3, 0.0);

    }
}
