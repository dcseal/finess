#include "dogdefs.h"

// This is a user-required routine.
//
//     Simple advection equation
//
// Input:
//
//    xpts( 1:numpts )          - The x-coordinates for a list of points
//    Q  ( 1:numpts, 1:meqn )   - The solution at each of these points
//    Aux( 1:numpts, 1:maux )   - The auxilary function at each point.
//
// Output:
//
//    flux( 1:numpts, 1:meqn )  - The flux function f(q) defined at each point
//
// See also: ...
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
