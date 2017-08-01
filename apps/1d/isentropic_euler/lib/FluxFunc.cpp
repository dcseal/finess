#include "dogdefs.h"
#include <math.h>
#include "IniParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Euler equations for gas dynamics
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{
  const int numpts=xpts.getsize();

  // Gas constant
  const double gamma = global_ini_params.get_gamma();

  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);      
    
      // Variables
      const double rho    = Q.get(i,1);
      const double u1     = Q.get(i,2)/rho;
      const double press=pow(rho,gamma); 
      // Flux function
      flux.set(i,1, rho*u1 );
      flux.set(i,2, rho*u1*u1 + press );
    }
}
