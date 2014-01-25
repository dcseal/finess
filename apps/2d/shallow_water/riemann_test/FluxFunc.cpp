#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     2D Shallow Water Equations
//
void FluxFunc(
        const dTensor2& xpts,
        const dTensor2& Q,
        const dTensor2& Aux,
        dTensor3& flux)
{
    double h,u1,u2;
    int i;
    int numpts=xpts.getsize(1);

    for (i=1; i<=numpts; i++)
    {
        // Variables
        h      = Q.get(i,1);
        u1     = Q.get(i,2)/h;
        u2     = Q.get(i,3)/h;

        // 1-component of flux function
        flux.set(i,1,1, h*u1 );
        flux.set(i,2,1, h*u1*u1 + 0.5*h*h );
        flux.set(i,3,1, h*u1*u2 );

        // 2-component of flux function
        flux.set(i,1,2, h*u2 );
        flux.set(i,2,2, h*u2*u1 );
        flux.set(i,3,2, h*u2*u2 + 0.5*h*h );
    }
}
