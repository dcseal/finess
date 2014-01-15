#include "dogdefs.h"

// This is a user-required routine.
//
// This routine defines the source-term function in a hyperbolic balance law.
//
// Input:
//
//    xpts( 1:numpts )              - The x-coordinates for a list of points
//    qvals  ( 1:numpts, 1:meqn )   - The solution at each of these points
//    auxvals( 1:numpts, 1:maux )   - The auxilary function at each point.
//
// Output:
//
//    Psi( 1:numpts, 1:meqn )       - The source term function psi(q,x) defined at each point
//
// See also: ...
//
void SourceTermFunc(const dTensor1& xpts, 
        const dTensor2& qvals, 
        const dTensor2& auxvals,
        dTensor2& Psi)
{

    const int numpts=xpts.getsize();
    const int meqn=qvals.getsize(2);
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        for (int me=1; me<=meqn; me++)
        {
            Psi.set(i,me, 0.0 );
        }
    }           

}
