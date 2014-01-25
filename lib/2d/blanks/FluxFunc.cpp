#include "tensors.h"

// * TEMPLATE *
//
// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
// Inputs:
//
//    xpts( 1:numpts, 1:ndim )  - ndim == 2 for 2D problems
//       Q( 1:numpts, 1:meqn )  - meqn == # of eqns in the system
//     Aux( 1:numpts, 1:maux )  - maux == # of aux variables used for the problem
// 
// Outputs:
//
//    flux( 1:numpts, 1:meqn, 1:ndim ) - ndim == 1: f, ndim == 2: g
//
// See also: blanks/*
void FluxFunc(const dTensor2& xpts, const dTensor2& Q,const dTensor2& Aux, dTensor2& flux)
{

    const int numpts=xpts.getsize(1);

    for(int i=1; i<=numpts; i++)
    for(int m=1; m<=meqn;   m++)
    {

        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        // Variables
        double qc = Q.get(i,m);

        // Flux function
        flux.set(i, m, 1, u*qc );
        flux.set(i, m, 2, v*qc );
    }

}
