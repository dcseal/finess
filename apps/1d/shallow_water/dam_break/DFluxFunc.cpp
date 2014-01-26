#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// Shallow Water Equations
//
void DFluxFunc(const dTensor1& xpts, 
	       const dTensor2& Q,
	       const dTensor2& Aux, 
	       dTensor3& Dflux)
{
  int i; //,j;
  int numpts=xpts.getsize();
//  double x;

  double h,u;

  for (i=1; i<=numpts; i++)
  {
//    x = xpts.get(i);
	
    h = Q.get(i,1);    // height
	u = Q.get(i,2)/h;  // velocity

    Dflux.set(i,1,1, 0.0e0);
	Dflux.set(i,1,2, 1.0e0);
	Dflux.set(i,2,1, -pow( u, 2) + h );
	Dflux.set(i,2,2, 2.0 * u);

  }

}
