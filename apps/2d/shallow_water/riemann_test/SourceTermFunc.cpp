#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc( const dTensor2& xpts, const dTensor2& qvals,
     const dTensor2& auxvals, dTensor2& fvals)
{
    int i,m;
    int numpts=xpts.getsize(1);
    int meqn=fvals.getsize(2);
    double x,y,h,b,bx,by;
	
    for (i=1; i<=numpts; i++)
    {
	x = xpts.get(i,1);
        y = xpts.get(i,2);

	h  = qvals.get(i,1);
	b  = auxvals.get(i,1);
	bx = auxvals.get(i,2);
	by = auxvals.get(i,3);

	fvals.set(i,1,  0.0  );
	fvals.set(i,2, -h*bx );
	fvals.set(i,3, -h*by );
    }
}
