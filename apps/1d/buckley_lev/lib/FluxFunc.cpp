#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Buckley-Leveret Equation
//
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{    

    const int numpts=xpts.getsize();

    // TODO - add this as a parameter supplied in parameters.ini
    const double M = (1./3.);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        double q  = Q.get(i,1);
        double q2 = Q.get(i,1)*Q.get(i,1);

        // Flux function
        flux.set(i,1, q2 / ( q2 + M*(1.0-q)*(1.0-q) ) );
    }

}
