#include <cmath>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const double gamma = global_ini_params.get_gamma();
    const int    numpts = xpts.getsize(1);

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);

        double rad = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        double rho,press,u1,u2,u3,energy;

        if(rad<0.3e0)
        {
            rho   =  global_ini_params.get_rhol();
            u1    =  global_ini_params.get_u1l();
            u2    =  global_ini_params.get_u2l();
            u3    =  global_ini_params.get_u3l();
            press =  global_ini_params.get_pl();
        }
        else
        {
            rho   =  global_ini_params.get_rhor();
            u1    =  global_ini_params.get_u1r();
            u2    =  global_ini_params.get_u2r();
            u3    =  global_ini_params.get_u3r();
            press =  global_ini_params.get_pr();
        }

        energy = press/(gamma-1.0e0)
            + 0.5e0*rho*(u1*u1+u2*u2+u3*u3);

        qvals.set(i,1, rho );       // density
        qvals.set(i,2, rho*u1 );    // 1-momentum
        qvals.set(i,3, rho*u2 );    // 2-momentum
        qvals.set(i,4, rho*u3 );    // 3-momentum
        qvals.set(i,5, energy );    // energy
    }

}
