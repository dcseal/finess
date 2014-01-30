#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// The expected format is Dflux.get(:, i, j) = \partial f_i, \partial q_j.
//
//     Simple advection equation, f'(q) = u
//     Burger's equation, f'(q) = q
//     Acoustics equation, f'(q) = [0 1; 1 0]
//
// Inputs:
//
//     xpts( 1:numpts )         - a list of x-values at various spatial points.
//        Q( 1:numpts, 1:meqn ) - a vector of conserved variables
//      Aux( 1:numpts, 1:maux ) - vector of auxilary values
//   
// Output:
//
//    Dflux( 1:numpts, 1:meqn, 1:meqn ) - f'(q) at each point.
//
// See also: FluxFunc and D2FluxFunc.
//
void DFluxFunc(const dTensor1& xpts, 
	       const dTensor2& Q,
	       const dTensor2& Aux,
	       dTensor3& Dflux)
{

    const int numpts=xpts.getsize();
    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        Dflux.set(i,1,1, 0. );
    }

}
