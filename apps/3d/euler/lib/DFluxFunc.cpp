#include <cmath>
#include "tensors.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// The expected format is Dflux.get(:,i,j,1:2) = 
//            (\partial f_i, \partial q_j, \partial g_i, \partial q_j )
//
//
void DFluxFunc(const dTensor2& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux,
        dTensor4& Dflux)
{
    using std::pow;
    const int numpts = xpts.getsize(1);
    const int meqn   = Q.getsize(2);

    // Gas constant
    const double gamma = global_ini_params.get_gamma();

    Dflux.setall(0.);
    for (int i=1; i<=numpts; i++)
    {

        // Not used: u3 ... 
        double q1 = Q.get(i,1);
        double q2 = Q.get(i,2);
        double q3 = Q.get(i,3);
        double q4 = Q.get(i,4);
        double q5 = Q.get(i,5);
        Dflux.set(i,1,1,1, 0);
        Dflux.set(i,1,2,1, 1);
        Dflux.set(i,1,3,1, 0);
        Dflux.set(i,1,4,1, 0);
        Dflux.set(i,1,5,1, 0);
        Dflux.set(i,2,1,1, (gamma*pow( q2, 2 ) + gamma*pow( q3, 2 ) + gamma*pow( q4, 2 ) - 3*pow( q2, 2 ) - pow( q3, 2 ) - pow( q4, 2 ))/(2*pow( q1, 2 )));
        Dflux.set(i,2,2,1, q2*(-gamma + 3)/q1);
        Dflux.set(i,2,3,1, q3*(-gamma + 1)/q1);
        Dflux.set(i,2,4,1, q4*(-gamma + 1)/q1);
        Dflux.set(i,2,5,1, gamma - 1);
        Dflux.set(i,3,1,1, -q2*q3/pow( q1, 2 ));
        Dflux.set(i,3,2,1, q3/q1);
        Dflux.set(i,3,3,1, q2/q1);
        Dflux.set(i,3,4,1, 0);
        Dflux.set(i,3,5,1, 0);
        Dflux.set(i,4,1,1, -q2*q4/pow( q1, 2 ));
        Dflux.set(i,4,2,1, q4/q1);
        Dflux.set(i,4,3,1, 0);
        Dflux.set(i,4,4,1, q2/q1);
        Dflux.set(i,4,5,1, 0);
        Dflux.set(i,5,1,1, q2*(-gamma*q1*q5 + gamma*pow( q2, 2 ) + gamma*pow( q3, 2 ) + gamma*pow( q4, 2 ) - pow( q2, 2 ) - pow( q3, 2 ) - pow( q4, 2 ))/pow( q1, 3 ));
        Dflux.set(i,5,2,1, (2*gamma*q1*q5 - 3*gamma*pow( q2, 2 ) - gamma*pow( q3, 2 ) - gamma*pow( q4, 2 ) + 3*pow( q2, 2 ) + pow( q3, 2 ) + pow( q4, 2 ))/(2*pow( q1, 2 )));
        Dflux.set(i,5,3,1, q2*q3*(-gamma + 1)/pow( q1, 2 ));
        Dflux.set(i,5,4,1, q2*q4*(-gamma + 1)/pow( q1, 2 ));
        Dflux.set(i,5,5,1, gamma*q2/q1);

        //Computing the Jacobian of the flux function, g'(q)
        Dflux.set(i,1,1,2, 0);
        Dflux.set(i,1,2,2, 0);
        Dflux.set(i,1,3,2, 1);
        Dflux.set(i,1,4,2, 0);
        Dflux.set(i,1,5,2, 0);
        Dflux.set(i,2,1,2, -q2*q3/pow( q1, 2 ));
        Dflux.set(i,2,2,2, q3/q1);
        Dflux.set(i,2,3,2, q2/q1);
        Dflux.set(i,2,4,2, 0);
        Dflux.set(i,2,5,2, 0);
        Dflux.set(i,3,1,2, (gamma*pow( q2, 2 ) + gamma*pow( q3, 2 ) + gamma*pow( q4, 2 ) - pow( q2, 2 ) - 3*pow( q3, 2 ) - pow( q4, 2 ))/(2*pow( q1, 2 )));
        Dflux.set(i,3,2,2, q2*(-gamma + 1)/q1);
        Dflux.set(i,3,3,2, q3*(-gamma + 3)/q1);
        Dflux.set(i,3,4,2, q4*(-gamma + 1)/q1);
        Dflux.set(i,3,5,2, gamma - 1);
        Dflux.set(i,4,1,2, -q3*q4/pow( q1, 2 ));
        Dflux.set(i,4,2,2, 0);
        Dflux.set(i,4,3,2, q4/q1);
        Dflux.set(i,4,4,2, q3/q1);
        Dflux.set(i,4,5,2, 0);
        Dflux.set(i,5,1,2, q3*(-gamma*q1*q5 + gamma*pow( q2, 2 ) + gamma*pow( q3, 2 ) + gamma*pow( q4, 2 ) - pow( q2, 2 ) - pow( q3, 2 ) - pow( q4, 2 ))/pow( q1, 3 ));
        Dflux.set(i,5,2,2, q2*q3*(-gamma + 1)/pow( q1, 2 ));
        Dflux.set(i,5,3,2, (2*gamma*q1*q5 - gamma*pow( q2, 2 ) - 3*gamma*pow( q3, 2 ) - gamma*pow( q4, 2 ) + pow( q2, 2 ) + 3*pow( q3, 2 ) + pow( q4, 2 ))/(2*pow( q1, 2 )));
        Dflux.set(i,5,4,2, q3*q4*(-gamma + 1)/pow( q1, 2 ));
        Dflux.set(i,5,5,2, gamma*q3/q1);

        //Computing the Jacobian of the flux function, h'(q)
        Dflux.set(i,1,1,3, 0);
        Dflux.set(i,1,2,3, 0);
        Dflux.set(i,1,3,3, 0);
        Dflux.set(i,1,4,3, 1);
        Dflux.set(i,1,5,3, 0);
        Dflux.set(i,2,1,3, -q2*q4/pow( q1, 2 ));
        Dflux.set(i,2,2,3, q4/q1);
        Dflux.set(i,2,3,3, 0);
        Dflux.set(i,2,4,3, q2/q1);
        Dflux.set(i,2,5,3, 0);
        Dflux.set(i,3,1,3, -q3*q4/pow( q1, 2 ));
        Dflux.set(i,3,2,3, 0);
        Dflux.set(i,3,3,3, q4/q1);
        Dflux.set(i,3,4,3, q3/q1);
        Dflux.set(i,3,5,3, 0);
        Dflux.set(i,4,1,3, (gamma*pow( q2, 2 ) + gamma*pow( q3, 2 ) + gamma*pow( q4, 2 ) - pow( q2, 2 ) - pow( q3, 2 ) - 3*pow( q4, 2 ))/(2*pow( q1, 2 )));
        Dflux.set(i,4,2,3, q2*(-gamma + 1)/q1);
        Dflux.set(i,4,3,3, q3*(-gamma + 1)/q1);
        Dflux.set(i,4,4,3, q4*(-gamma + 3)/q1);
        Dflux.set(i,4,5,3, gamma - 1);
        Dflux.set(i,5,1,3, q4*(-gamma*q1*q5 + gamma*pow( q2, 2 ) + gamma*pow( q3, 2 ) + gamma*pow( q4, 2 ) - pow( q2, 2 ) - pow( q3, 2 ) - pow( q4, 2 ))/pow( q1, 3 ));
        Dflux.set(i,5,2,3, q2*q4*(-gamma + 1)/pow( q1, 2 ));
        Dflux.set(i,5,3,3, q3*q4*(-gamma + 1)/pow( q1, 2 ));
        Dflux.set(i,5,4,3, (2*gamma*q1*q5 - gamma*pow( q2, 2 ) - gamma*pow( q3, 2 ) - 3*gamma*pow( q4, 2 ) + pow( q2, 2 ) + pow( q3, 2 ) + 3*pow( q4, 2 ))/(2*pow( q1, 2 )));
        Dflux.set(i,5,5,3, gamma*q4/q1);


    }

}
