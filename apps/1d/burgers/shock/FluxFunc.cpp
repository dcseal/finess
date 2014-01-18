#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Burgers equation
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{    
    int i;
    int numpts=xpts.getsize();
    double x;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);

        // Flux function
        flux.set(i,1, 0.5*pow(Q.get(i,1),2) );
    }

}
