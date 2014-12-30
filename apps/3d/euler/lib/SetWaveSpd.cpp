#include "dogdefs.h"
#include "dog_math.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
void SetWaveSpd(const dTensor1& nvec, 
		const dTensor1& xedge,
		const dTensor1& Ql, 
		const dTensor1& Qr,
		const dTensor1& Auxl, 
		const dTensor1& Auxr,
		double& s1,
		double& s2)
{   
  // Gas constant
  double const gamma   = global_ini_params.get_gamma();

  // Left states
  double const rhol    = Ql.get(1);
  double const ul      = (nvec.get(1)*Ql.get(2)+nvec.get(2)*Ql.get(3)+nvec.get(3)*Ql.get(4))/rhol;
  double const energyl = Ql.get(5);
  double const pressl  = (gamma-1.0e0)*(energyl-0.5e0*(pow(Ql.get(2),2)+pow(Ql.get(3),2)+pow(Ql.get(4),2))/rhol);

  // Right states
  double const rhor    = Qr.get(1);
  double const ur      = (nvec.get(1)*Qr.get(2)+nvec.get(2)*Qr.get(3)+nvec.get(3)*Qr.get(4))/rhor;
  double const energyr = Qr.get(5);
  double const pressr  = (gamma-1.0e0)*(energyr-0.5e0*(pow(Qr.get(2),2)+pow(Qr.get(3),2)+pow(Qr.get(4),2))/rhor);

  // Average states
  double const rho    = 0.5e0*(rhol+rhor);
  double const u      = 0.5e0*(ul+ur);
  double const press  = 0.5e0*(pressl+pressr);

  // Sound speeds
  double const cl = sqrt(fabs(gamma*pressl/rhol));
  double const cr = sqrt(fabs(gamma*pressr/rhor));
  double const c  = sqrt(fabs(gamma*press/rho));
  
  // Minimum speed
  s1 = Min(ul-cl, u-c);
  
  // Maximum speed
  s2 = Max(ur+cr, u+c);
  
}
