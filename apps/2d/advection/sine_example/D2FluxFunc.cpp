#include "dogdefs.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
// In 2d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     Simple Advection Equation
//
void D2FluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor4& D2flux)
{
    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);

        // 2nd Derivative of Flux function
        D2flux.set(i, 1, 1, 1, 0.0e0 );
        D2flux.set(i, 1, 1, 2, 0.0e0 );

    }    
}
