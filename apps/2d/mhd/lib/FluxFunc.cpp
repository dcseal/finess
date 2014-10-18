#include "tensors.h"

#include "IniParams.h"
// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     2D Ideal MHD equations
//
void FluxFunc(const dTensor2& xpts,
	      const dTensor2& Q,
	      const dTensor2& Aux,
	      dTensor3& flux)
{
  const int numpts=xpts.getsize(1);

  // Gas constant
  double const gamma = global_ini_params.get_gamma();

  for (int i=1; i<=numpts; i++){
    // Variables
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
      
    // 1-component of flux function
    flux.set(i,1,1, rho*u1 );
    flux.set(i,2,1, rho*u1*u1 + press +  Bm - B1*B1 );
    flux.set(i,3,1, rho*u1*u2 - B1*B2 );
    flux.set(i,4,1, rho*u1*u3 - B1*B3 );
    flux.set(i,5,1, u1*(energy + press + Bm) - B1*Bu );
    flux.set(i,6,1, 0.0 );
    flux.set(i,7,1, u1*B2 - u2*B1 );
    flux.set(i,8,1, u1*B3 - u3*B1 );

    // 2-component of flux function
    flux.set(i,1,2, rho*u2);
    flux.set(i,2,2, rho*u1*u2 - B2*B1);
    flux.set(i,3,2, rho*u2*u2 + press + Bm - B2*B2);
    flux.set(i,4,2, rho*u2*u3 - B2*B3);
    flux.set(i,5,2, u2*(energy + press + Bm)  -  B2*Bu);
    flux.set(i,6,2, u2*B1 - u1*B2);
    flux.set(i,7,2, 0.0);
    flux.set(i,8,2, u2*B3 - u3*B2);
  }
}
