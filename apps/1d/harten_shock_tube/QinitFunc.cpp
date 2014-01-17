#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor1& xpts, 
	       dTensor2& qvals)
{
    const int numpts=xpts.getsize();

    // Gas constant
    const double gamma = eulerParams.gamma;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        double rho,m1,E;

        if(x<0.5e0)
        {
            rho   = 0.445;
            m1    = 0.3111;
            E     = 8.928;
        }
        else
        {
            rho   = 0.5;
            m1    = 0.0;
            E     = 1.4275;
        }

        qvals.set(i,1, rho );       // density
        qvals.set(i,2, m1  );       // 1-momentum
        qvals.set(i,3, 0.0 );    // 2-momentum
        qvals.set(i,4, 0.0 );    // 3-momentum
        qvals.set(i,5, E   );    // energy		
    }
}
