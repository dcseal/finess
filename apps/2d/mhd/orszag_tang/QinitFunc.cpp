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
    using std::pow;
    using std::sin;

    const int numpts=xpts.getsize(1);

    const double gamma_gas = global_ini_params.get_gamma();
    const double gamma1 = global_ini_params.get_gamma() - 1.0;
    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i,1);
      double y = xpts.get(i,2);
      double rho,u1,u2,u3,press,B1,B2,B3,energy;


      rho   =  pow(gamma_gas, 2);
      u1    =  -sin(y);
      u2    =  sin(x);
      u3    =  0.0;
      press =  gamma_gas;
      B1    =  -sin(y);
      B2    =  sin(2.0*x);
      B3    =  0.0;



      energy = press/gamma1 
	+ 0.5*rho*(u1*u1 + u2*u2 + u3*u3)
	+ 0.5*(B1*B1 + B2*B2 + B3*B3);

      qvals.set(i,1, rho );     // density
      qvals.set(i,2, rho*u1);  // 1-momentum
      qvals.set(i,3, rho*u2 );  // 2-momentum
      qvals.set(i,4, rho*u3 );  // 3-momentum
      qvals.set(i,5, energy );  // energy
      qvals.set(i,6,  B1 );      // B1
      qvals.set(i,7,  B2 );      // B2
      qvals.set(i,8, B3 );      // B3      
    }
}
