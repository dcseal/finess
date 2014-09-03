#include "tensors.h"

// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
		    const dTensor2& qvals, 
		    const dTensor2& auxvals,
                    dTensor2& fvals)
{
    int i;
    int numpts=xpts.getsize();
    double x;
    
    for (i=1; i<=numpts; i++)
    {
	x = xpts.get(i);

	fvals.set(i,1,  0.0 );
    }
}
