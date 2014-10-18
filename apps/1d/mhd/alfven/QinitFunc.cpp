#include <cmath>
#include <fstream>
#include "tensors.h"
#include "IniParams.h"
using namespace std;

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
void QinitFunc(const dTensor1& xpts, 
        dTensor2& qvals)
{
    const int numpts=xpts.getsize();
    const double gamma1 = global_ini_params.get_gamma()-1.0;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        double rho,u1,u2,u3,press,B1,B2,B3,energy;

        rho   =  1.0;
        u1    =  0.0;
        u2    =  0.1 * sin(2.0*M_PI * x);
        u3    =  0.1 * cos(2.0*M_PI * x);
        press =  0.1;
        B1    =  1.0;
        B2    =  0.1 * sin(2.0*M_PI * x);
        B3    =  0.1 * cos(2.0*M_PI * x);
        
        energy = press/gamma1 
            + 0.5*rho*(u1*u1 + u2*u2 + u3*u3)
            + 0.5*(B1*B1 + B2*B2 + B3*B3);

        qvals.set(i,1, rho );     // density
        qvals.set(i,2, rho*u1 );  // 1-momentum
        qvals.set(i,3, rho*u2 );  // 2-momentum
        qvals.set(i,4, rho*u3 );  // 3-momentum
        qvals.set(i,5, energy );  // energy
        qvals.set(i,6, B1 );      // B1
        qvals.set(i,7, B2 );      // B2
        qvals.set(i,8, B3 );      // B3      
    }
}
