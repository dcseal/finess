#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Compressible Euler equations
//
void FluxFunc(const dTensor2& xpts,
        const dTensor2& Q, 
        const dTensor2& Aux, 
        dTensor3& flux)
{

    const int   numpts = xpts.getsize(1);
    const double gamma = global_ini_params.get_gamma();

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);

        // Variables
        const double   rho  = Q.get(i,1);
        const double    u1  = Q.get(i,2)/rho;
        const double    u2  = Q.get(i,3)/rho;
        const double    u3  = Q.get(i,4)/rho;      
        const double energy = Q.get(i,5);
        const double press  = (gamma-1.0)*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

        // Flux function: x-direction 
        flux.set(i,1,1, rho*u1             );
        flux.set(i,2,1, rho*u1*u1 + press  );
        flux.set(i,3,1, rho*u1*u2          );
        flux.set(i,4,1, rho*u1*u3          );
        flux.set(i,5,1, u1*(energy+press)  );

        // Flux function: y-direction 
        flux.set(i,1,2, rho*u2             );
        flux.set(i,2,2, rho*u2*u1          );
        flux.set(i,3,2, rho*u2*u2 + press  );
        flux.set(i,4,2, rho*u2*u3          );
        flux.set(i,5,2, u2*(energy+press)  );

        // Flux function: z-direction 
        flux.set(i,1,3, rho*u3             );
        flux.set(i,2,3, rho*u3*u1          );
        flux.set(i,3,3, rho*u3*u2          );
        flux.set(i,4,3, rho*u3*u3 + press  );
        flux.set(i,5,3, u3*(energy+press)  );
    }

}
