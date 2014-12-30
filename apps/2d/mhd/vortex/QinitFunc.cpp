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
    using std::exp;

    const int numpts=xpts.getsize(1);

    const double gamma_gas = global_ini_params.get_gamma();
    const double gamma1 = global_ini_params.get_gamma() - 1.0;
    const double kappa = global_ini_params.get_kappa();
    const double mu = global_ini_params.get_mu();

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
      const double rho_m = 1.0, press_m = 1.0, u1_m = 1.0, u2_m = 1.0, u3_m = 0.0;
      double x = xpts.get(i,1);
      double y = xpts.get(i,2);
      double rho,u1,u2,u3,press,B1,B2,B3,energy;

      double r2 = x*x + y*y;
      double exphalfsomething = exp(0.5 * (1.0 - r2));
      double expsomething = exp(1.0 - r2);
      
      rho = rho_m;
      press = press_m + 1.0/(8.0*M_PI*M_PI) * expsomething
                         *((1-r2)*mu*mu - kappa*kappa);
      u1 = u1_m + kappa/(2.0*M_PI) * exphalfsomething * (-y);
      u2 = u2_m + kappa/(2.0*M_PI) * exphalfsomething * x;
      u3 = u3_m;

      B1 = mu/(2.0*M_PI) * exphalfsomething * (-y);
      B2 = mu/(2.0*M_PI) * exphalfsomething * x;
      B3 = 0.0;

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
