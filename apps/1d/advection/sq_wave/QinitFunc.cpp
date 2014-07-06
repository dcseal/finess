#include <cmath>
#include "dogdefs.h"

// This is a user-required routine.
//
// This routine defines the initial conditions for the problem.
//
// Input:
//
//    xpts( 1:numpts, 1:ndim )
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )
//
// See also: ...
//
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize();
    const double width = 0.5;

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        if( fabs( x ) > width )
        { qvals.set(i,1, 0.0e0 ); }
        else
        { qvals.set(i, 1,  1.0 ); }

    }

}
