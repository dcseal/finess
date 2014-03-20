#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize();

    // Gas constant
    const double gamma = eulerParams.gamma;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        double rho = 1.0;
        double u1  = 0.0;
        double u2  = 0.0;
        double u3  = 0.0;
        double press = 1.0e2;

        if (x<0.1)
        {  press = 1.0e3;  }
        else if (x<0.9)
        {  press = 1.0e-2;  }        

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1+u2*u2+u3*u3);

        qvals.set(i,1, rho );       // density
        qvals.set(i,2, rho*u1 );    // 1-momentum
        qvals.set(i,3, rho*u2 );    // 2-momentum
        qvals.set(i,4, rho*u3 );    // 3-momentum
        qvals.set(i,5, energy );    // energy 
    }
}
