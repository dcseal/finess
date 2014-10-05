#include "tensors.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Ideal MHD equations
//
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{
  const int numpts=xpts.getsize();

  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);

      // Gas constant
      const double gamma = global_ini_params.get_gamma();

      // Variables
      const double rho    = Q.get(i,1);
      const double u1     = Q.get(i,2)/rho;
      const double u2     = Q.get(i,3)/rho;
      const double u3     = Q.get(i,4)/rho;
      const double energy = Q.get(i,5);
      const double B1     = Q.get(i,6);
      const double B2     = Q.get(i,7);
      const double B3     = Q.get(i,8);
      const double press  = (gamma-1.0e0)*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3)
			      -0.5*(B1*B1 + B2*B2 + B3*B3));
      const double Bm     = 0.5*(B1*B1 + B2*B2 + B3*B3);
      const double Bu     = u1*B1 + u2*B2 + u3*B3;  
      
      // Flux function
      flux.set(i,1, rho*u1 );
      flux.set(i,2, rho*u1*u1 + press +  Bm - B1*B1 );
      flux.set(i,3, rho*u1*u2 - B1*B2 );
      flux.set(i,4, rho*u1*u3 - B1*B3 );
      flux.set(i,5, u1*(energy + press + Bm) - B1*Bu );
      flux.set(i,6, 0.0 );
      flux.set(i,7, u1*B2 - u2*B1 );
      flux.set(i,8, u1*B3 - u3*B1 );
    }
}
