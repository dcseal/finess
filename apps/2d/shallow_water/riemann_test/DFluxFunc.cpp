#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
//     Shallow Water equations
//
void DFluxFunc(const dTensor2& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor4& Dflux)
{

    const int numpts = xpts.getsize(1);
    const int meqn   = Q.getsize(2);

//  Dflux.setall(0.);
    for (int i=1; i<=numpts; i++)
    {

        double q1 = Q.get(i,1);
        double q2 = Q.get(i,2);
        double q3 = Q.get(i,3);

        //////////////////////
        // Jacobian, f'(q): //
        //////////////////////

        Dflux.set(i,1,1,1, 0.0);
        Dflux.set(i,1,2,1, 1.0);
        Dflux.set(i,1,3,1, 0.0);
        Dflux.set(i,2,1,1, (pow(q1,3) - pow(q2,2))/pow(q1,2));
        Dflux.set(i,2,2,1, 2.0*q2/q1);
        Dflux.set(i,2,3,1, 0.0);
        Dflux.set(i,3,1,1, -q2*q3/pow(q1,2));
        Dflux.set(i,3,2,1, q3/q1);
        Dflux.set(i,3,3,1, q2/q1);
 
        //////////////////////
        // Jacobian, g'(q): //
        //////////////////////
        Dflux.set(i,1,1,2, 0.0);
        Dflux.set(i,1,2,2, 0.0);
        Dflux.set(i,1,3,2, 1.0);
        Dflux.set(i,2,1,2, -q2*q3/pow(q1,2));
        Dflux.set(i,2,2,2, q3/q1);
        Dflux.set(i,2,3,2, q2/q1);
        Dflux.set(i,3,1,2, (pow(q1,3) - pow(q3,2))/pow(q1,2));
        Dflux.set(i,3,2,2, 0.0);
        Dflux.set(i,3,3,2, 2.0*q3/q1);

    }

}
