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
//    xpts ( 1:numpts, 1:ndim )  - The x,y, and z-coordinates for a list of points
//   qvals ( 1:numpts, 1:meqn )
//  auxvals( 1:numpts, 1:maux )
//
// Output:
//
//    source( 1:numpts, 1:meqn ) - The source function for each equation.
//
// See also: AuxFunc.
void SourceTermFunc(const dTensor2& xpts, 
                    const dTensor2& qvals, 
                    const dTensor2& auxvals,
                    dTensor2& source )
{

    const int numpts = xpts.getsize(1);
    const int meqn   = source.getsize(2);

    for(int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);

        for(int m=1; m<=meqn; m++)
        {
            source.set(i,m, 0.0e0 );
        }
    }
}
