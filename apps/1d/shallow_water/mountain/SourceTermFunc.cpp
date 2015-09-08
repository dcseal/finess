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
		
	h = qvals.get(i,1);
	u = qvals.get(i,2)/h;
	
	bot  = auxvals.get(i,1);
	dbdx = auxvals.get(i,2);
	
	fvals.set(i,1,  0.0 );
	fvals.set(i,2, -h*dbdx );
    }
}
