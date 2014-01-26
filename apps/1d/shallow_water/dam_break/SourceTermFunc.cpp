#include "tensors.h"

// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
		    const dTensor2& qvals, 
		    const dTensor2& auxvals,
                    dTensor2& fvals)
{
    int i,me;
    int numpts=xpts.getsize();
    int meqn=qvals.getsize(2);
    double x,h,u,bot,dbdx;
    
    for (i=1; i<=numpts; i++)
    {
	x = xpts.get(i);
	
	for (me=1; me<=meqn; me++)
	{  fvals.set(i,me, 0.0 );  }
    }
}
