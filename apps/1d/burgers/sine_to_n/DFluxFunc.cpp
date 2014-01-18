#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
// 
// The expected format is Dflux.get(:, i, j) = \partial f_i, \partial q_j.
//
//     Burger's equation, f'(q) = q
//
void DFluxFunc(const dTensor1& xpts, 
	       const dTensor2& Q,
	       const dTensor2& Aux, 
	       dTensor3& Dflux)
{
  int i,j;
  int numpts=xpts.getsize();
  double x;

  for (i=1; i<=numpts; i++)
  {
    x = xpts.get(i);

	//Dflux.set(i,1,1, Aux.get(i,1)  );  // Simple Advection Equation
	Dflux.set(i,1,1, Q.get(i,1)  ); // Burger's equation
  }    
}
