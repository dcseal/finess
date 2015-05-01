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
    //const double epsilon = 10.0828;
    const double epsilon = 5.0;
    const double t0 = 0.0;

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1)-t0*1.0;
        double y = xpts.get(i,2)-t0*1.0;
        if (x<-5.0) x = x + 10.0;
        if (y<-5.0) y = y + 10.0;

        double rho   =  1.0;
        double u1    =  1.0;
        double u2    =  1.0;
        double u3    =  0.0;
        double press =  1.0;
        double tpr   = press/rho;
        const double entpy = press/pow(rho,gamma);

        double r2 = x*x+y*y;
        double du = -epsilon/2.e0/pi*exp((1.e0-r2)/2.e0)*y;
        double dv =  epsilon/2.e0/pi*exp((1.e0-r2)/2.e0)*x;
        double dtpr = (1.e0-gamma)/gamma * epsilon*epsilon/8.e0/pi/pi*exp(1.e0-r2);

        u1 = u1 + du;
        u2 = u2 + dv;

        tpr = tpr + dtpr;
        press=pow(tpr,gamma/(gamma-1.e0)) / pow(entpy,1.e0/(gamma-1.e0));
        rho = press/tpr;

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
