#include "tensors.h"
#include "EulerParams.h"
#include "constants.h"
#include <cmath>

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

        double rho = 1.0 + 0.5*sin(2.0*pi*x);
        double u1  = 1.0;
        double u2  = 0.0;
        double u3  = 0.0;
        double press = 1.0;

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1+u2*u2+u3*u3);

        qvals.set(i,1, rho );       // density
        qvals.set(i,2, rho*u1 );    // 1-momentum
        qvals.set(i,3, rho*u2 );    // 2-momentum
        qvals.set(i,4, rho*u3 );    // 3-momentum
        qvals.set(i,5, energy );    // energy 
    }
}
