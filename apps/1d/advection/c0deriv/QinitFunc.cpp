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

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        qvals.set(i, 1,  0.0 );
        if( 0.25 < x & x < 0.4 )
        { qvals.set(i,1, (x-0.25)/0.075); }
        else if( 0.4 <= x & x < 0.6 )
        { qvals.set(i, 1, 2.0 ); }
        else if( 0.6 <= x & x < 0.75 )
        { qvals.set(i,1, (0.75-x)/0.075); }

    }

}
