#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Ideal MHD equations
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
                const dTensor1& Auxr,
		double& s1,double& s2)
{   
  // Gas constant
  const double gamma = global_ini_params.get_gamma();
    
  // Left states
  const double rhol = Ql.get(1);
  const double u1l  = Ql.get(2)/rhol;
  const double u2l  = Ql.get(3)/rhol;
  const double u3l  = Ql.get(4)/rhol;
  const double energyl = Ql.get(5);
  const double B1l  = Ql.get(6);
  const double B2l  = Ql.get(7);
  const double B3l  = Ql.get(8);
  const double Bm2l = B1l*B1l + B2l*B2l + B3l*B3l;
  const double um2l = u1l*u1l + u2l*u2l + u3l*u3l;
  const double pressl = (gamma-1.0e0)*(energyl  -0.5*rhol*um2l - 0.5*Bm2l);
  
  // Right states
  const double rhor = Qr.get(1);
  const double u1r  = Qr.get(2)/rhor;
  const double u2r  = Qr.get(3)/rhor;
  const double u3r  = Qr.get(4)/rhor;
  const double energyr = Qr.get(5);
  const double B1r  = Qr.get(6);
  const double B2r  = Qr.get(7);
  const double B3r  = Qr.get(8);
  const double Bm2r = B1r*B1r + B2r*B2r + B3r*B3r;
  const double um2r = u1r*u1r + u2r*u2r + u3r*u3r;
  const double pressr = (gamma-1.0e0)*(energyr - 0.5*rhor*um2r - 0.5*Bm2r);
  
  // Average states
  const double rho    = 0.5e0*(rhol+rhor);
  const double u1     = 0.5e0*(u1l+u1r);
  const double u2     = 0.5e0*(u2l+u2r);
  const double u3     = 0.5e0*(u3l+u3r);
  const double press  = 0.5e0*(pressl+pressr);
  const double B1     = 0.5e0*(B1l+B1r);
  const double B2     = 0.5e0*(B2l+B2r);
  const double B3     = 0.5e0*(B3l+B3r);
  const double Bm2    = B1*B1 + B2*B2 + B3*B3;
  
  // Sound speed squared
  const double a2l = (gamma*pressl/rhol);
  const double a2r = sqrt(fabs(gamma*pressr/rhor));
  const double a2  = sqrt(fabs(gamma*press/rho));
  
  // Fast magnetosonic speeds
  const double cf_l = sqrt(0.5*(a2l + Bm2l/rhol + 
		   sqrt(fabs(pow(a2l+Bm2l/rhol,2) - 4.0*a2l*B1l*B1l/rhol))));
  const double cf_r = sqrt(0.5*(a2r + Bm2r/rhor + 
		   sqrt(fabs(pow(a2r+Bm2r/rhor,2) - 4.0*a2r*B1r*B1r/rhor))));
  const double cf   = sqrt(0.5*(a2 + Bm2/rho + 
		   sqrt(fabs(pow(a2+Bm2/rho,2) - 4.0*a2*B1*B1/rho))));
  
  // Minimum speed
  s1 = Min(u1l - cf_l, u1 - cf);
  
  // Maximum speed
  s2 = Max(u1r + cf_r, u1 + cf);
    
}
