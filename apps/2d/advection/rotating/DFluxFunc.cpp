#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// The expected format is Dflux.get(:,i,j,1:2) = 
//            (\partial f_i, \partial q_j, \partial g_i, \partial q_j )
//
//     Simple advection equation, f'(q) = u1, g'(q) = u2.
//     Burger's equation, f'(q) = q, g'(q) = q
//
void DFluxFunc(const dTensor2& xpts, 
	       const dTensor2& Q,
	       const dTensor2& Aux,
	       dTensor4& Dflux)
{

    const int numpts=xpts.getsize(1);
    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        // Variables 
        double qc = Q.get(i,1);
        double u  = Aux.get(i,1);
        double v  = Aux.get(i,2);

        Dflux.set(i,1,1,1, u );  // First component of flux func
        Dflux.set(i,1,1,2, v );  // Second component of flux func
    }

}
