#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
// 
// The expected format is Dflux.get(:, i, j) = \partial f_i, \partial q_j.
//
void DFluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor3& Dflux)
{

    const int numpts=xpts.getsize();

    // TODO - add this as a parameter in parameters.ini ...
    const double M = (1./3.);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        double q  = Q.get(i,1);
        double q2 = Q.get(i,1)*Q.get(i,1);

        double den = q2 + M*(1.0-q)*(1.0-q);
        Dflux.set(i,1,1, 2.0*M*q*(1.0-q)/ (den*den) );
    }    
}
