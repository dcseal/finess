#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
// In 2d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     Simple Advection Equation
//
void D2FluxFunc(const dTensor2& xpts, 
		const dTensor2& Q,
		const dTensor2& Aux,
		dTensor5& D2flux)
{

    const int numpts = xpts.getsize(1);
    const int meqn   = Q.getsize(2);

    D2flux.setall(0.);
    for(int i=1; i<=numpts; i++){
        double q = Q.get(i, 1);
        D2flux.set(i, 1, 1, 1, 1, 1.0);
        D2flux.set(i, 1, 1, 1, 2, 1.0);
        D2flux.set(i, 1, 1, 1, 3, 1.0);
    }

}
