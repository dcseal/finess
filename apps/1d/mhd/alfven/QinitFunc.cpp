#include <fstream>
#include "tensors.h"
#include "MHDParams.h"
#include "InitialParams.h"
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
    const double gamma1 = mhdParams.gamma-1.0;

    const double rhol = initialParams.rhol;
    const double u1l  = initialParams.unl;
    const double u2l  = initialParams.utl;
    const double u3l  = initialParams.u3l;
    const double pl   = initialParams.pl;
    const double B1l  = initialParams.Bnl;
    const double B2l  = initialParams.Btl;
    const double B3l  = initialParams.B3l;

    const double rhor = initialParams.rhor;
    const double u1r  = initialParams.unr;
    const double u2r  = initialParams.utr;
    const double u3r  = initialParams.u3r;
    const double pr   = initialParams.pr;
    const double B1r  = initialParams.Bnr;
    const double B2r  = initialParams.Btr;
    const double B3r  = initialParams.B3r;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        double rho,u1,u2,u3,press,B1,B2,B3,energy;

        if(x<0.0e0)
        {
            rho   =  rhol;
            u1    =  u1l;
            u2    =  u2l;
            u3    =  u3l;
            press =  pl;
            B1    =  B1l;
            B2    =  B2l;
            B3    =  B3l;
        }
        else
        {
            rho   =  rhor;
            u1    =  u1r;
            u2    =  u2r;
            u3    =  u3r;
            press =  pr;
            B1    =  B1r;
            B2    =  B2r;
            B3    =  B3r;
        }

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
