#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
// In 2d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     Shallow Water equations
//
void D2FluxFunc(const dTensor2& xpts, 
		const dTensor2& Q,
		const dTensor2& Aux,
		dTensor5& D2flux)
{

    const int numpts = xpts.getsize(1);
    const int meqn   = Q.getsize(2);

    D2flux.setall(0.);

    for (int i=1; i<=numpts; i++)
    {

        double q1 = Q.get(i,1);
        double q2 = Q.get(i,2);
        double q3 = Q.get(i,3);

        //////////////////////
        // Hessian, f''(q): //
        //////////////////////

//      D2flux.set(i,1,1,1,1, 0);
//      D2flux.set(i,1,1,2,1, 0);
//      D2flux.set(i,1,1,3,1, 0);
//      D2flux.set(i,1,2,1,1, 0);
//      D2flux.set(i,1,2,2,1, 0);
//      D2flux.set(i,1,2,3,1, 0);
//      D2flux.set(i,1,3,1,1, 0);
//      D2flux.set(i,1,3,2,1, 0);
//      D2flux.set(i,1,3,3,1, 0);
         
        D2flux.set(i,2,1,1,1, (pow(q1,3) + 2*pow(q2,2))/pow(q1,3));
        D2flux.set(i,2,1,2,1, -2*q2/pow(q1,2));
        D2flux.set(i,2,1,3,1, 0);
        D2flux.set(i,2,2,1,1, -2*q2/pow(q1,2));
        D2flux.set(i,2,2,2,1, 2/q1);
        D2flux.set(i,2,2,3,1, 0);
        D2flux.set(i,2,3,1,1, 0);
        D2flux.set(i,2,3,2,1, 0);
        D2flux.set(i,2,3,3,1, 0);
         
        D2flux.set(i,3,1,1,1, 2*q2*q3/pow(q1,3));
        D2flux.set(i,3,1,2,1, -q3/pow(q1,2));
        D2flux.set(i,3,1,3,1, -q2/pow(q1,2));
        D2flux.set(i,3,2,1,1, -q3/pow(q1,2));
        D2flux.set(i,3,2,2,1, 0);
        D2flux.set(i,3,2,3,1, 1/q1);
        D2flux.set(i,3,3,1,1, -q2/pow(q1,2));
        D2flux.set(i,3,3,2,1, 1/q1);
        D2flux.set(i,3,3,3,1, 0);
         
        //////////////////////
        // Hessian, g''(q): //
        //////////////////////

//      D2flux.set(i,1,1,1,2, 0);
//      D2flux.set(i,1,1,2,2, 0);
//      D2flux.set(i,1,1,3,2, 0);
//      D2flux.set(i,1,2,1,2, 0);
//      D2flux.set(i,1,2,2,2, 0);
//      D2flux.set(i,1,2,3,2, 0);
//      D2flux.set(i,1,3,1,2, 0);
//      D2flux.set(i,1,3,2,2, 0);
//      D2flux.set(i,1,3,3,2, 0);
         
        D2flux.set(i,2,1,1,2, 2*q2*q3/pow(q1,3));
        D2flux.set(i,2,1,2,2, -q3/pow(q1,2));
        D2flux.set(i,2,1,3,2, -q2/pow(q1,2));
        D2flux.set(i,2,2,1,2, -q3/pow(q1,2));
        D2flux.set(i,2,2,2,2, 0);
        D2flux.set(i,2,2,3,2, 1/q1);
        D2flux.set(i,2,3,1,2, -q2/pow(q1,2));
        D2flux.set(i,2,3,2,2, 1/q1);
        D2flux.set(i,2,3,3,2, 0);
         
        D2flux.set(i,3,1,1,2, (pow(q1,3) + 2*pow(q3,2))/pow(q1,3));
        D2flux.set(i,3,1,2,2, 0);
        D2flux.set(i,3,1,3,2, -2*q3/pow(q1,2));
        D2flux.set(i,3,2,1,2, 0);
        D2flux.set(i,3,2,2,2, 0);
        D2flux.set(i,3,2,3,2, 0);
        D2flux.set(i,3,3,1,2, -2*q3/pow(q1,2));
        D2flux.set(i,3,3,2,2, 0);
        D2flux.set(i,3,3,3,2, 2/q1);

    }

}
