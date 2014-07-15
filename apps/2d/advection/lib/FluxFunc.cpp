#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Simple advection equation
//
void FluxFunc(const dTensor2& xpts, const dTensor2& Q, 
        const dTensor2& Aux, dTensor3& flux)
{

    const int numpts=xpts.getsize(1);
    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        // Variables
        double qc = Q.get(i,1);
        double u  = Aux.get(i,1);
        double v  = Aux.get(i,2);

        // Flux function
        flux.set(i, 1, 1, u*qc );
        flux.set(i, 1, 2, v*qc );
    }

}
