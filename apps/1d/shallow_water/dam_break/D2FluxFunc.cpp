#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
//     Shallow Water Equations
//
void D2FluxFunc(const dTensor1& xpts, 
		const dTensor2& Q,
		const dTensor2& Aux, 
		dTensor4& D2flux)
{
  int i;
  int numpts=xpts.getsize();
//  double x;
  double h,u;

  for (i=1; i<=numpts; i++)
  {
//    x = xpts.get(i);
	h = Q.get(i,1);
	u = Q.get(i,2)/h;
	
    D2flux.set(i, 1, 1, 1, 0.0);
    D2flux.set(i, 1, 1, 2, 0.0);
    D2flux.set(i, 1, 2, 1, 0.0);
    D2flux.set(i, 1, 2, 2, 0.0);
    D2flux.set(i, 2, 1, 1, 1.0 + 2.0*u*u/h );
    D2flux.set(i, 2, 1, 2, - 2.0*u/h);
    D2flux.set(i, 2, 2, 1, - 2.0*u/h);
    D2flux.set(i, 2, 2, 2, 2.0 * pow(h,-1) );

  }

}
