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
    const double width = 2.0*0.25;

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        double s = x;//-0.5;
        
        if ( fabs(s) > width/2.0 )
        { qvals.set(i,1, 0.0e0 ); }
        else
        { qvals.set(i, 1,  pow(cos(pi*s/width),6) ); }
        qvals.set(i,2,0.0e0);

    }
  
}
