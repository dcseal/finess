#include <cmath>
#include "tensors.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// Euler Equations
//
void DFluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor3& Dflux)
{

    const int numpts=xpts.getsize();

    // Gas constant
    const double gamma = global_ini_params.get_gamma();

    Dflux.setall(0.);
    for (int i=1; i<=numpts; i++)
    {

        // Not used: u2 and u3 ... 
        const double q1 = Q.get(i,1);
        const double q2 = Q.get(i,2);
        const double q3 = Q.get(i,3);
        const double q4 = Q.get(i,4);
        const double q5 = Q.get(i,5);

        const double rho    = Q.get(i,1);
        const double u1     = Q.get(i,2)/rho;
        const double u2     = Q.get(i,3)/rho;
        const double u3     = Q.get(i,4)/rho;
        const double en     = Q.get(i,5);
        const double press  = (gamma-1.0e0)*(en-0.5e0*rho*(u1*u1+u2*u2+u3*u3));
 
        // Computing the Jacobian of the flux function, f'(q)
        
        // pd{ f_1 }{ q_j }
        Dflux.set(i,1,1, 0);
        Dflux.set(i,1,2, 1.00000000000000);
        Dflux.set(i,1,3, 0);
        Dflux.set(i,1,4, 0);
        Dflux.set(i,1,5, 0);

        // pd{ f_2 }{ q_j }
        Dflux.set(i,2,1, 0.5*(gamma*pow( q2, 2 ) + gamma*pow( q3, 2 ) + gamma*pow( q4, 2 ) - 3.0*pow( q2, 2 ) - pow( q3, 2 ) - pow( q4, 2 ))/pow( q1, 2 ));
        Dflux.set(i,2,2, q2*(-gamma + 3.0)/q1);
        Dflux.set(i,2,3, q3*(-gamma + 1.0)/q1);
        Dflux.set(i,2,4, q4*(-gamma + 1.0)/q1);
        Dflux.set(i,2,5, 0);

        // pd{ f_3 }{ q_j }
        Dflux.set(i,3,1, -q2*q3/pow( q1, 2 ));
        Dflux.set(i,3,2, q3/q1);
        Dflux.set(i,3,3, q2/q1);
        Dflux.set(i,3,4, 0);
        Dflux.set(i,3,5, 0);

        // pd{ f_4 }{ q_j }
        Dflux.set(i,4,1, -q2*q4/pow( q1, 2 ));
        Dflux.set(i,4,2, q4/q1);
        Dflux.set(i,4,3, 0);
        Dflux.set(i,4,4, q2/q1);
        Dflux.set(i,4,5, 0);

        // pd{ f_5 }{ q_j }
        Dflux.set(i,5,1, q2*(-en*gamma*q1 + en*q1 + gamma*pow( q2, 2 ) + gamma*pow( q3, 2 ) + gamma*pow( q4, 2 ) - q1*q5 - pow( q2, 2 ) - pow( q3, 2 ) - pow( q4, 2 ))/pow( q1, 3 ));
        Dflux.set(i,5,2, 0.5*(2.0*en*gamma*q1 - 2.0*en*q1 - 3.0*gamma*pow( q2, 2 ) - gamma*pow( q3, 2 ) - gamma*pow( q4, 2 ) + 2.0*q1*q5 + 3.0*pow( q2, 2 ) + pow( q3, 2 ) + pow( q4, 2 ))/pow( q1, 2 ));
        Dflux.set(i,5,3, q2*q3*(-gamma + 1.0)/pow( q1, 2 ));
        Dflux.set(i,5,4, q2*q4*(-gamma + 1.0)/pow( q1, 2 ));
        Dflux.set(i,5,5, q2/q1);

    }

}
