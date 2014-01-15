#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Simple advection equation
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{    
    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {

        // Flux function f(q) in q_t + f(q)_x = psi.
        double x = xpts.get(i);
        flux.set(i,1, Q.get(i,1) * Aux.get(i,1) );

    }

}
