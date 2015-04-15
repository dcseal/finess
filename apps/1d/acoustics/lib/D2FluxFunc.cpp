#include "dogdefs.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
//     Simple Advection Equation
//
void D2FluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor4& D2flux)
{
    const int numpts=xpts.getsize();
    const int meqn=2;
    for (int i=1; i<=numpts; i++)
    for (int m1=1; m1<=meqn; m1++)
    for (int m2=1; m2<=meqn; m2++)
    for (int m3=1; m3<=meqn; m3++)
    {
        D2flux.set(i, m1, m2, m3, 0.0);
    }


}
