#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);
    double h,u1,u2,b;

    // Loop over grid points
    for(int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

	if( sqrt( pow(x-0.4,2)+pow(y-0.4,2) ) < 0.35)
	{
	    h   =  1.0;
	    u1  =  0.0;
	    u2  =  0.0;
	}
	else
	{
	    h   =  0.1;
	    u1  =  0.0;
	    u2  =  0.0;
	} 

        qvals.set(i, 1, h    );
        qvals.set(i, 2, h*u1 );
        qvals.set(i, 3, h*u2 );
    }
}
