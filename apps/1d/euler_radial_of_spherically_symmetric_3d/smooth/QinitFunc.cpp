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
    using std::exp;
    using std::pow;
    const int numpts=xpts.getsize();

    // Gas constant
    const double gamma = global_ini_params.get_gamma();

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        double rho = 1.0 + 0.1*exp(-30.0 * pow(x-1.0, 2));
        double u1  = 0.0;
//        double press = 1.0;

//        double energy = press/(gamma-1.0e0) 
//            + 0.5e0*rho*(u1*u1+u2*u2+u3*u3);
        double energy = rho;
        qvals.set(i,1, rho );       // density
        qvals.set(i,2, rho*u1 );    // 1-momentum
        qvals.set(i,3, energy );    // energy 
    }
}
