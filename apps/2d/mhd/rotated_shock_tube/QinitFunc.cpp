#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "MHDParams.h"
#include "InitialParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

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

    const double angle = initialParams.angle;

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i,1);
      double y = xpts.get(i,2);
      double rho,u1,u2,u3,press,B1,B2,B3,energy;

      double xi = x*cos(angle) + y*sin(angle);
      if(xi<0)
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
      qvals.set(i,2, rho*( u1*cos(angle) + u2*sin(angle)) );  // 1-momentum
      qvals.set(i,3, rho*(-u1*sin(angle) + u2*cos(angle)) );  // 2-momentum
      qvals.set(i,4, rho*u3 );  // 3-momentum
      qvals.set(i,5, energy );  // energy
      qvals.set(i,6,  B1*cos(angle) + B2*sin(angle) );      // B1
      qvals.set(i,7, -B1*sin(angle) + B2*cos(angle) );      // B2
      qvals.set(i,8, B3 );      // B3      
    }
}
