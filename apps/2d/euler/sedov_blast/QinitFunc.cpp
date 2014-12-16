#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// For Sedov blast problem, this is supposed to be a delta function in energy
// at x=y=0.  To approximate this, we simply set a single spike at the lower
// left most corner of the domain.
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

    // Gas consant
    const double gamma = global_ini_params.get_gamma();
    const double dx    = global_ini_params.get_dx();
    const double dy    = global_ini_params.get_dy();

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);
        double rho,u1,u2,u3,energy;

        rho = 1.0;
        u1  = 0.0;
        u2  = 0.0;
        u3  = 0.0;
        if(x < dx - 1.e-4 && y < dy - 1.e-4 )
        { energy =  0.244816/dx/dy; }
        else
        { energy =  1.0e-12; } 

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
