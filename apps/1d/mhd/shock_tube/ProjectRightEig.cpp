#include <cmath>
#include "tensors.h"
#include "MHDParams.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(const dTensor1& Aux_ave, 
		     const dTensor1& Q_ave, 
		     const dTensor2& Wvals,
		     dTensor2& Qvals)
{    
  const int meqn = Qvals.getsize(1);
  const int kmax = Qvals.getsize(2)+1;
  const int ixy = 1;  
  dTensor2 ru(8,8),Rmat(8,8);
  
  // Directional information
  int mu1,mu2,mb1,mb2;
  if (ixy==1)
    {
      mu1 = 2;
      mu2 = 3;
      mb1 = 6;
      mb2 = 7;
    }
  
  // Average states
  const double gamma  = mhdParams.gamma;
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
  
  const double a2 = fabs(gamma*press/rho);
  const double a = sqrt(a2);
  const double d = a2 + (B1*B1 + B2*B2 + B3*B3)/rho;
  
  const double ca = sqrt(fabs(B1*B1/rho));
  const double cf = sqrt(0.5*( d + sqrt(d*d - 4.0*a2*B1*B1/rho)));
  const double cs = sqrt(0.5*( d - sqrt(d*d - 4.0*a2*B1*B1/rho)));
  
  double beta1,beta2,beta3,alphaf,alphas;

  if (B1<0.0)
    {  beta1 = -1.0; }
  else
    {  beta1 =  1.0;  }

  if ( (B2*B2 + B3*B3) <= 1.0e-12 ) 
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
  
  //     # 1 - right eigenvector
  ru.set(1,1,    rho*alphaf );
  ru.set(mu1,1, -alphaf*cf );
  ru.set(mu2,1,  alphas*beta2*cs*beta1 );
  ru.set(4,1,    alphas*beta3*cs*beta1 );
  ru.set(5,1,    alphaf*gamma*press );
  ru.set(mb1,1,  0.0 );
  ru.set(mb2,1,  alphas*beta2*a*rhosq );
  ru.set(8,1,    alphas*beta3*a*rhosq );
  
  //     # 2 - right eigenvector
  ru.set(1,2,    0.0 );
  ru.set(mu1,2,  0.0 );
  ru.set(mu2,2, -beta3/twosq );
  ru.set(4,2,    beta2/twosq );
  ru.set(5,2,    0.0 );
  ru.set(mb1,2,  0.0 );
  ru.set(mb2,2, -beta3*rhosq/twosq );
  ru.set(8,2,    beta2*rhosq/twosq );
  
  //     # 3 - right eigenvector
  ru.set(1,3,    rho*alphas );
  ru.set(mu1,3, -alphas*cs );
  ru.set(mu2,3, -alphaf*beta2*cf*beta1 );
  ru.set(4,3,   -alphaf*beta3*cf*beta1 );
  ru.set(5,3,    alphas*gamma*press );
  ru.set(mb1,3,  0.0 );
  ru.set(mb2,3, -alphaf*beta2*a*rhosq );
  ru.set(8,3,   -alphaf*beta3*a*rhosq );
  
  //     # 4 - right eigenvector
  ru.set(1,4,    1.0 );
  ru.set(mu1,4,  0.0 );
  ru.set(mu2,4,  0.0 );
  ru.set(4,4,    0.0 );
  ru.set(5,4,    0.0 );
  ru.set(mb1,4,  0.0 );
  ru.set(mb2,4,  0.0 );
  ru.set(8,4,    0.0 );
  
  //     # 5 - right eigenvector
  ru.set(1,5,    0.0 );
  ru.set(mu1,5,  0.0 );
  ru.set(mu2,5,  0.0 );
  ru.set(4,5,    0.0 );
  ru.set(5,5,    0.0 );
  ru.set(mb1,5,  1.0 );
  ru.set(mb2,5,  0.0 );
  ru.set(8,5,    0.0 );
  
  //     # 6 - right eigenvector
  ru.set(1,6,    rho*alphas );
  ru.set(mu1,6,  alphas*cs );
  ru.set(mu2,6,  alphaf*beta2*cf*beta1 );
  ru.set(4,6,    alphaf*beta3*cf*beta1 );
  ru.set(5,6,    alphas*gamma*press );
  ru.set(mb1,6,  0.0 );
  ru.set(mb2,6, -alphaf*beta2*a*rhosq );
  ru.set(8,6,   -alphaf*beta3*a*rhosq );
  
  //     # 7 - right eigenvector
  ru.set(1,7,    0.0 );
  ru.set(mu1,7,  0.0 );
  ru.set(mu2,7, -beta3/twosq );
  ru.set(4,7,    beta2/twosq );
  ru.set(5,7,    0.0 );
  ru.set(mb1,7,  0.0 );
  ru.set(mb2,7,  beta3*rhosq/twosq );
  ru.set(8,7,   -beta2*rhosq/twosq );

  //     # 8 - right eigenvector
  ru.set(1,8,    rho*alphaf );
  ru.set(mu1,8,  alphaf*cf );
  ru.set(mu2,8, -alphas*beta2*cs*beta1 );
  ru.set(4,8,   -alphas*beta3*cs*beta1 );
  ru.set(5,8,    alphaf*gamma*press );
  ru.set(mb1,8,  0.0 );
  ru.set(mb2,8,  alphas*beta2*a*rhosq );
  ru.set(8,8,    alphas*beta3*a*rhosq );
  
  //     # ---------------------------------------------------
  //     # CONSERVATIVE E-VECTORS
  //     # ---------------------------------------------------
  for (int m=1; m<=8; m++)
    {
      Rmat.set(1,m,    ru.get(1,m) );
      Rmat.set(mu1,m,  ru.get(mu1,m)*rho + ru.get(1,m)*u1 );
      Rmat.set(mu2,m,  ru.get(mu2,m)*rho + ru.get(1,m)*u2 );
      Rmat.set(4,m,    ru.get(4,m)*rho   + ru.get(1,m)*u3 );
      Rmat.set(5,m,    um2*ru.get(1,m) + rho*u1*ru.get(mu1,m) 
	       + rho*u2*ru.get(mu2,m) + rho*u3*ru.get(4,m) + ru.get(5,m)/gm1 
	       + B1*ru.get(mb1,m) + B2*ru.get(mb2,m) + B3*ru.get(8,m) );
      Rmat.set(mb1,m,  ru.get(mb1,m) );
      Rmat.set(mb2,m,  ru.get(mb2,m) );
      Rmat.set(8,m,    ru.get(8,m) );
    }
  
  // Project onto right eigenvectors
  for (int k=1; k<=(kmax-1); k++)
    for (int m=1; m<=8; m++)
      {
	Qvals.set(m,k, 0.0 );
	
	for (int ell=1; ell<=8; ell++)
	  {
	    double tmp = Qvals.get(m,k);
	    Qvals.set(m,k, tmp + Rmat.get(m,ell)*Wvals.get(ell,k) );
	  }
      }
}
