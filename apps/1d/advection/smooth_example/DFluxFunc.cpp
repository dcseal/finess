#include "dogdefs.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
//     Simple advection equation, f'(q) = u
//
void DFluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor3& Dflux)
{
    const int numpts = xpts.getsize();
    const int meqn   = Q.getsize(2);

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        Dflux.set(i,1,1, Aux.get(i,1)  );
    }    
}
