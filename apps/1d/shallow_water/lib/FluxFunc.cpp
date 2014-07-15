#include "tensors.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Shallow water equations
//
void FluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor2& flux)
{

    const int numpts=xpts.getsize();

    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        // Variables
        double h = Q.get(i,1);
        double u = Q.get(i,2)/h;

        // Flux function
        flux.set(i,1, h*u );
        flux.set(i,2, h*u*u + 0.5*h*h );
    }    
}
