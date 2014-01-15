#include "dogdefs.h"

// This is a user-supplied routine that sets the source term
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
