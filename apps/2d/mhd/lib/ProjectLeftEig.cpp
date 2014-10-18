#include <cmath>
#include "dogdefs.h"

#include "IniParams.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(int ixy, 
		    const dTensor1& Aux_ave,
		    const dTensor1& Q_ave, 
		    const dTensor2& Qvals, 
		    dTensor2& Wvals)
{    
  
  const int meqn = Qvals.getsize(1);
  const int kmax = Qvals.getsize(2)+1;
  dTensor2 ru(8,8),Lmat(8,8);  

  // Directional information
  int mu1,mu2,mb1,mb2;
  if (ixy==1) {
    mu1 = 2;
    mu2 = 3;
    mb1 = 6;
    mb2 = 7;
  }
  else{
    mu1 = 3;
    mu2 = 2;
    mb1 = 7;
    mb2 = 6;
  }
  
  // Average states
  const double gamma  = global_ini_params.get_gamma();
  const double rho    = Q_ave.get(1);
  const double u1     = Q_ave.get(2)/rho;
  const double u2     = Q_ave.get(3)/rho;
  const double u3     = Q_ave.get(4)/rho;
  const double energy = Q_ave.get(5);
  const double B1     = Q_ave.get(6);
  const double B2     = Q_ave.get(7);
  const double B3     = Q_ave.get(8);
  const double um2    = 0.5*(u1*u1 + u2*u2 + u3*u3);
  const double Bm2    = 0.5*(B1*B1 + B2*B2 + B3*B3);
  const double press  = (gamma-1.0)*(energy - rho*um2 - Bm2);
    
  // Some useful quantities
  const double rhosq = sqrt(rho);
  const double twosq = sqrt(2.0);
  
  const double gm1 = gamma-1.0;
  
  const double a2 = (gamma*press/rho);
  const double a = sqrt(a2);
  const double d = a2 + (B1*B1 + B2*B2 + B3*B3)/rho;
  
  const double ca = sqrt(B1*B1/rho);
  const double cf = sqrt(0.5*( d + sqrt(d*d - 4.0*a2*B1*B1/rho)));
  const double cs = sqrt(0.5*( d - sqrt(d*d - 4.0*a2*B1*B1/rho)));

  double beta1,beta2,beta3,alphaf,alphas;

  if (B1<0.0)
    {  beta1 = -1.0; }
  else
    {  beta1 =  1.0;  }
  
  if ( (B2*B2 + B3*B3)<=1.0e-12 )
    {
      beta2 = 1.0/twosq;
      beta3 = 1.0/twosq;
    }
  else
    {
      beta2 = B2/sqrt(B2*B2 + B3*B3);
      beta3 = B3/sqrt(B2*B2 + B3*B3);
    }
  
  if ( (B2*B2 + B3*B3)<=1.0e-12 &&
       fabs(a2 - B1*B1/rho)<=1.0e-12 )
    {
      alphaf = sin(atan(1.0)/2.0);
      alphas = cos(atan(1.0)/2.0);
    }
  else
    {
      alphaf = sqrt(fabs(a2 - cs*cs))/sqrt(fabs(cf*cf-cs*cs));
      alphas = sqrt(fabs(cf*cf - a2))/sqrt(fabs(cf*cf-cs*cs));
    }
  
  // ---------------------------------------------------
  // PRIMITIVE E-VECTORS
  // ---------------------------------------------------
  
  //     # 1 - left eigenvector
  ru.set(1,1,    0.0 );
  ru.set(1,mu1, -alphaf*cf/(2.0*a2) );
  ru.set(1,mu2,  alphas*cs*beta2*beta1/(2.0*a2) );
  ru.set(1,4,    alphas*cs*beta3*beta1/(2.0*a2) );
  ru.set(1,5,    alphaf/(2.0*rho*a2) );
  ru.set(1,mb1,  0.0 );
  ru.set(1,mb2,  alphas*beta2/(2.0*rhosq*a) );
  ru.set(1,8,    alphas*beta3/(2.0*rhosq*a) );
  
  //     # 2 - left eigenvector
  ru.set(2,1,    0.0 );
  ru.set(2,mu1,  0.0 );
  ru.set(2,mu2, -beta3/twosq );
  ru.set(2,4,    beta2/twosq );
  ru.set(2,5,    0.0 );
  ru.set(2,mb1,  0.0 );
  ru.set(2,mb2, -beta3/(twosq*rhosq) );
  ru.set(2,8,    beta2/(twosq*rhosq) );
    
  //     # 3 - left eigenvector
  ru.set(3,1,    0.0 );
  ru.set(3,mu1, -alphas*cs/(2.0*a2) );
  ru.set(3,mu2, -alphaf*cf*beta2*beta1/(2.0*a2) );
  ru.set(3,4,   -alphaf*cf*beta3*beta1/(2.0*a2) );
  ru.set(3,5,    alphas/(2.0*rho*a2) );
  ru.set(3,mb1,  0.0 );
  ru.set(3,mb2, -alphaf*beta2/(2.0*rhosq*a) );
  ru.set(3,8,   -alphaf*beta3/(2.0*rhosq*a) );
  
  //     # 4 - left eigenvector
  ru.set(4,1,    1.0 );
  ru.set(4,mu1,  0.0 );
  ru.set(4,mu2,  0.0 );
  ru.set(4,4,    0.0 );
  ru.set(4,5,   -1.0/a2 );
  ru.set(4,mb1,  0.0 );
  ru.set(4,mb2,  0.0 );
  ru.set(4,8,    0.0 );
  
  //     # 5 - left eigenvector
  ru.set(5,1,    0.0 );
  ru.set(5,mu1,  0.0 );
  ru.set(5,mu2,  0.0 );
  ru.set(5,4,    0.0 );
  ru.set(5,5,    0.0 );
  ru.set(5,mb1,  1.0 );
  ru.set(5,mb2,  0.0 );
  ru.set(5,8,    0.0 );
  
  //     # 6 - left eigenvector
  ru.set(6,1,    0.0 );
  ru.set(6,mu1,  alphas*cs/(2.0*a2) );
  ru.set(6,mu2,  alphaf*cf*beta2*beta1/(2.0*a2) );
  ru.set(6,4,    alphaf*cf*beta3*beta1/(2.0*a2) );
  ru.set(6,5,    alphas/(2.0*rho*a2) );
  ru.set(6,mb1,  0.0 );
  ru.set(6,mb2, -alphaf*beta2/(2.0*rhosq*a) );
  ru.set(6,8,   -alphaf*beta3/(2.0*rhosq*a) );
  
  //     # 7 - left eigenvector
  ru.set(7,1,    0.0 );
  ru.set(7,mu1,  0.0 );
  ru.set(7,mu2, -beta3/twosq );
  ru.set(7,4,    beta2/twosq );
  ru.set(7,5,    0.0 );
  ru.set(7,mb1,  0.0 );
  ru.set(7,mb2,  beta3/(twosq*rhosq) );
  ru.set(7,8,   -beta2/(twosq*rhosq) );
  
  //     # 8 - left eigenvector
  ru.set(8,1,    0.0 );
  ru.set(8,mu1,  alphaf*cf/(2.0*a2) );
  ru.set(8,mu2, -alphas*cs*beta2*beta1/(2.0*a2) );
  ru.set(8,4,   -alphas*cs*beta3*beta1/(2.0*a2) );
  ru.set(8,5,    alphaf/(2.0*rho*a2) );
  ru.set(8,mb1,  0.0 );
  ru.set(8,mb2,  alphas*beta2/(2.0*rhosq*a) );
  ru.set(8,8,    alphas*beta3/(2.0*rhosq*a) );

  // ---------------------------------------------------
  // CONSERVATIVE E-VECTORS
  // ---------------------------------------------------
  for (int m=1; m<=8; m++)
    {
      Lmat.set(m,1,     ru.get(m,1)-ru.get(m,mu1)*u1/rho
	       -ru.get(m,mu2)*u2/rho
	       -ru.get(m,4)*u3/rho+ru.get(m,5)*um2*gm1 );
      Lmat.set(m,mu1,   ru.get(m,mu1)/rho-ru.get(m,5)*u1*gm1 );
      Lmat.set(m,mu2,   ru.get(m,mu2)/rho-ru.get(m,5)*u2*gm1 );
      Lmat.set(m,4,     ru.get(m,4)/rho  -ru.get(m,5)*u3*gm1 );
      Lmat.set(m,5,     ru.get(m,5)*gm1 );
      Lmat.set(m,mb1,  -ru.get(m,5)*B1*gm1+ru.get(m,mb1) );
      Lmat.set(m,mb2,  -ru.get(m,5)*B2*gm1+ru.get(m,mb2) );
      Lmat.set(m,8,    -ru.get(m,5)*B3*gm1+ru.get(m,8) );
    }
  
  // Project onto left eigenvectors
  for (int k=1; k<=(kmax-1); k++)    
    for (int m=1; m<=8; m++)
      {
	Wvals.set(m,k, 0.0 );
        
	for (int ell=1; ell<=8; ell++)
	  {
	    double tmp = Wvals.get(m,k);
	    Wvals.set(m,k, tmp + Lmat.get(m,ell)*Qvals.get(ell,k) );
	  }
      } 

}
