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
    double h,u;
    
    // Initial conditions
    for (i=1; i<=numpts; i++)
    {
	x = xpts.get(i);
	
	if (x<0.5)
	{
	    h = 3.0;
	    u = 0.0;
	}
	else
	{
	    h = 1.0;
	    u = 0.0;
	}

        qvals.set(i,1, h );    // depth
        qvals.set(i,2, h*u );  // momentum	
    }
}
