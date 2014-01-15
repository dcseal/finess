#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor1& xpts, 
        dTensor2& auxvals)
{

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);

        // advection speed u, in q_t + (u q)_x = 0
        // This could be modified to contain spatial dependence.
        auxvals.set(i,1, 1.0e0 );
    }

}
