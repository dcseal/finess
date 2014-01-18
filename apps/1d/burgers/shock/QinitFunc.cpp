#include "tensors.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor1& xpts, 
	       dTensor2& qvals)
{
    int i;
    int numpts=xpts.getsize();
    double x;
    
    for (i=1; i<=numpts; i++)
    {
	x = xpts.get(i);
	if( x < 0.1 )
	{
		qvals.set(i,1,1.0);
	}else
	{
		qvals.set(i,1,0.0);
	}
	
    }
}
