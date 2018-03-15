#include <cmath>
#include "tensors.h"
#include "constants.h"
#include "IniParams.h"

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
    const double gamma = global_ini_params.get_gamma();

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        double eps=0.2;

        /*
        double rho = 0.955 + 5.0*eps*eps*(1.0-cos(2.0*pi*x));
        double sign = 1.0;
        if(x<0.0){sign=-sign;}
        double u1  = -10.0*eps*eps*sign*sqrt(1.4)*(1.0-cos(2.0*pi*x));
        */
        double rho = 1.0 + 0.2*sin(4.0*pi*x);
        double u1  = 1.0;

        double u2  = 0.0;
        double u3  = 0.0;
        double press = pow(rho,gamma);

        qvals.set(i,1, rho );       // density
        qvals.set(i,2, rho*u1 );    // 1-momentum
    }
}
