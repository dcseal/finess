#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// The expected format is Dflux.get(:,i,j,1:2) = 
//            (\partial f_i, \partial q_j, \partial g_i, \partial q_j )
//
//     Simple advection equation, f'(q) = u1, g'(q) = u2, h'(q) = u3.
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
        double z = xpts.get(i,3);

        // Variables 
        double qc = Q.get(i,1);
        double u1 = Aux.get(i,1);
        double u2 = Aux.get(i,2);
        double u3 = Aux.get(i,3);

        Dflux.set(i,1,1,1, u1 );  // First component of flux func
        Dflux.set(i,1,1,2, u2 );  // Second component of flux func
        Dflux.set(i,1,1,3, u3 );  // Third component of flux func
    }

}
