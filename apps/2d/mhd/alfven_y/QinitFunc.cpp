#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"

#include "IniParams.h"
// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

    const double gamma1 = global_ini_params.get_gamma() - 1.0;

//    const double angle = global_ini_params.get_angle();
    const double angle = M_PI / 2.0;
    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i,1);
      double y = xpts.get(i,2);
      double rho,u1,u2,u3,press,B1,B2,B3,energy;

//      double xi = x*cos(angle) + y*sin(angle);
      double xi = y;
      rho   =  1.0;
      u1    =  0.0;
      u2    =  0.1 * sin(2.0*M_PI * xi);
      u3    =  0.1 * cos(2.0*M_PI * xi);
      press =  0.1;
      B1    =  1.0;
      B2    =  0.1 * sin(2.0*M_PI * xi);
      B3    =  0.1 * cos(2.0*M_PI * xi);



      energy = press/gamma1 
	+ 0.5*rho*(u1*u1 + u2*u2 + u3*u3)
	+ 0.5*(B1*B1 + B2*B2 + B3*B3);

      qvals.set(i,1, rho );     // density
      qvals.set(i,2, rho*( u1*cos(angle) - u2*sin(angle)) );  // 1-momentum
      qvals.set(i,3, rho*( u1*sin(angle) + u2*cos(angle)) );  // 2-momentum
      qvals.set(i,4, rho*u3 );  // 3-momentum
      qvals.set(i,5, energy );  // energy
      qvals.set(i,6,  B1*cos(angle) - B2*sin(angle) );      // B1
      qvals.set(i,7,  B1*sin(angle) + B2*cos(angle) );      // B2
      qvals.set(i,8, B3 );      // B3      
    }
}
