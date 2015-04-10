#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
// In 1d this is a 3-tensor
//
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
// Buckley-Leverett Equation
//
void D2FluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor4& D2flux)
{
    int i,m1,m2,m3;
    int numpts=xpts.getsize();

    double x;
    double h,u;

    const double M = 1./3.;
    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);

        const double q2 = Q.get(i,1)*Q.get(i,1);
        const double q3 = q2*Q.get(i,1);

        const double qm1_sqd = (1.-Q.get(i,1))*(1.-Q.get(i,1));

        // 2nd Derivative of Flux function
        D2flux.set(i, 1, 1, 1, 
            2.*M*( -3.*q2 + M - 3.*M*q2 + 2.*q3*( 1.+M) ) * 
                pow( q2+M*qm1_sqd, -3 ) );

    }    
}
