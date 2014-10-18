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

    // Gas consant
    const double gamma = global_ini_params.get_gamma();
    const double gm1   = global_ini_params.get_gamma() - 1.0;
    const double gp1   = global_ini_params.get_gamma() + 1.0;

    const double Mach  = 5.09;
    const double M2    = Mach*Mach;

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double rho,u1,u2,u3,press;

        rho   =  1.4;
        u1    =  0.0;
        u2    =  0.0;
        u3    =  0.0;
        press =  1.0;

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        // Correct with the Mach number for incoming shock
        const double vs = Mach*sqrt( gamma*press/rho );
        if( x < 0.5 )
        {
            press = press*(2.0*gamma*M2- gm1)/gp1;
            rho   = rho*gp1*M2 / (gm1*M2 + 2.0);
            u1    = vs*(1.0 - 1.4 / rho );
            energy = press/(gamma-1.0e0) 
                + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);
        }

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }

}
