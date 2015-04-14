#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Burgers equation: f(q) = 0.5 q^2.
//
void FluxFunc(const dTensor1& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor2& flux)
{    
    const int numpts = xpts.getsize();
    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        // Flux function
        flux.set(i,1, 0.5*Q.get(i,1)*Q.get(i,1) );
    }

}
