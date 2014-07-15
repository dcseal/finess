#include "tensors.h"

// *REQUIRED* 
//
// This is a user-required routine that defines the initial conditions for the
// problem.
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts, 1:ndim )   - The x,y, and z-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);

        {  qvals.set(i,1, 0.0 ); }

    }
}
