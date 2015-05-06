#include <cmath>
#include "dogdefs.h"
#include "gas05.h"
#include "FiveMomentParams.h"

FiveMomentParams fiveMomentParams;
using namespace FiveMomentComponentID;

void FiveMomentFluxFunc(
    int n_offset,
    const dTensor2& Q,
    dTensor3& flux)
{
    // Parameters
    double gamma = fiveMomentParams.get_gamma();

    const int n_rho = n_offset + _rho;
    const int n_M1  = n_offset + _M1 ;
    const int n_M2  = n_offset + _M2 ;
    const int n_M3  = n_offset + _M3 ;
    const int n_N   = n_offset + _N  ;
        
    for(int j=1;j<=Q.getsize(1);j++)
    {
    
        // Variables
        const double& rho    = Q.get(j,n_rho);
        const double& M1     = Q.get(j,n_M1 );
        const double& M2     = Q.get(j,n_M2 );
        const double& M3     = Q.get(j,n_M3 );
        const double& energy = Q.get(j,n_N  );
        
        const double u1     = M1/rho;
        const double u2     = M2/rho;
        const double u3     = M3/rho;
    
        const double press  = (gamma-1.0e0)
            *(energy - 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3));
    
        // 1-component of flux function
        flux.set(j,n_rho,1,   M1 );
        flux.set(j,n_M1 ,1,   u1*M1 + press );
        flux.set(j,n_M2 ,1,   u1*M2 );
        flux.set(j,n_M3 ,1,   u1*M3 );
        flux.set(j,n_N  ,1,   u1*(energy+press) ); 
        
        // 2-component of flux function
        flux.set(j,n_rho,2,   M2 );
        flux.set(j,n_M1 ,2,   u2*M1 );
        flux.set(j,n_M2 ,2,   u2*M2 + press );
        flux.set(j,n_M3 ,2,   u2*M3 );
        flux.set(j,n_N  ,2,   u2*(energy+press) ); 
    }
}

void FiveMomentFluxFunc1(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
{
    // Parameters
    double gamma = fiveMomentParams.get_gamma();

    const int n_rho = n_offset + _rho;
    const int n_M1  = n_offset + _M1 ;
    const int n_M2  = n_offset + _M2 ;
    const int n_M3  = n_offset + _M3 ;
    const int n_N   = n_offset + _N  ;
        
    for(int j=1;j<=Q.getsize(1);j++)
    {
    
        // Variables
        const double& rho    = Q.get(j,n_rho);
        const double& M1     = Q.get(j,n_M1 );
        const double& M2     = Q.get(j,n_M2 );
        const double& M3     = Q.get(j,n_M3 );
        const double& energy = Q.get(j,n_N  );
        
        const double u1     = M1/rho;
        const double u2     = M2/rho;
        const double u3     = M3/rho;
    
        const double press  = (gamma-1.0e0)
            *(energy - 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3));
    
        // 1-component of flux function
        flux.set(j,n_rho, M1 );
        flux.set(j,n_M1 , u1*M1 + press );
        flux.set(j,n_M2 , u1*M2 );
        flux.set(j,n_M3 , u1*M3 );
        flux.set(j,n_N  , u1*(energy+press) ); 
    }
}

void FiveMomentFluxFunc2(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
{
    // Parameters
    double gamma = fiveMomentParams.get_gamma();

    const int n_rho = n_offset + _rho;
    const int n_M1  = n_offset + _M1 ;
    const int n_M2  = n_offset + _M2 ;
    const int n_M3  = n_offset + _M3 ;
    const int n_N   = n_offset + _N  ;
        
    for(int j=1;j<=Q.getsize(1);j++)
    {
    
        // Variables
        const double& rho    = Q.get(j,n_rho);
        const double& M1     = Q.get(j,n_M1 );
        const double& M2     = Q.get(j,n_M2 );
        const double& M3     = Q.get(j,n_M3 );
        const double& energy = Q.get(j,n_N  );
        
        const double u1     = M1/rho;
        const double u2     = M2/rho;
        const double u3     = M3/rho;
    
        const double press  = (gamma-1.0e0)
            *(energy - 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3));
    
        // 2-component of flux function
        flux.set(j,n_rho, M2 );
        flux.set(j,n_M1 , u2*M1 );
        flux.set(j,n_M2 , u2*M2 + press );
        flux.set(j,n_M3 , u2*M3 );
        flux.set(j,n_N  , u2*(energy+press) ); 
    }
}

void ProjectLeftEig_FiveMoment( int ixy, int n_offset,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{
    int mu,mv,mw;

    // Directional information
    if (ixy==1)
    {
	mu = 2;
	mv = 3;
	mw = 4;
    }
    else
    {
        assert(ixy==2);
	mu = 3;
	mv = 2;
	mw = 4;
    }

    const int n_rho = n_offset + _rho;
    const int nu    = n_offset + mu;
    const int nv    = n_offset + mv;
    const int nw    = n_offset + mw;
    const int n_N   = n_offset + _N;

    // Average states
    const double& gamma  = fiveMomentParams.get_gamma();
    const double& rho    = Q_ave.get(n_rho);
    const double  u1     = Q_ave.get(nu)/rho;
    const double  u2     = Q_ave.get(nv)/rho;
    const double  u3     = Q_ave.get(nw)/rho;
    const double& energy = Q_ave.get(n_N);
    const double  umag2  = (u1*u1 + u2*u2 + u3*u3);
    const double  press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
    const double  c      = sqrt(gamma*press/rho);
    const double  H      = (energy+press)/rho; 
  
    // Project onto left eigenvectors
    for (int k=1; k<=W.getsize(2); k++)
    {
        W.set(n_rho,k, ((umag2/2.0-H-u1*c)*Q.get(nu,k) 
                        + (umag2/2.0*(c-u1)+H*u1)*Q.get(n_rho,k) 
                        + c*(Q.get(n_N,k)-u2*Q.get(nv,k)
                        - u3*Q.get(nw,k)))/(c*(2.0*H-umag2)) );

        W.set(nu,k, 2.0*((H-umag2)*Q.get(n_rho,k) + u1*Q.get(nu,k) 
                        + u2*Q.get(nv,k) + u3*Q.get(nw,k) 
                        - Q.get(n_N,k))/(2.0*H-umag2) );

        W.set(nv,k, Q.get(nv,k)-u2*Q.get(n_rho,k) );

        W.set(nw,k, Q.get(nw,k)-u3*Q.get(n_rho,k) );

        W.set(n_N, k, ((H-umag2/2.0-u1*c)*Q.get(nu,k) 
                        + (umag2/2.0*(c+u1)-H*u1)*Q.get(n_rho,k)
                        + c*(Q.get(n_N,k)-u2*Q.get(nv,k) 
                        - u3*Q.get(nw,k)))/(c*(2.0*H-umag2)) );
    }
}

void ProjectRightEig_FiveMoment(int ixy, int n_offset,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{
    int mu,mv,mw;

    // Directional information
    if (ixy==1)
    {
	mu = 2;
	mv = 3;
	mw = 4;
    }
    else
    {
        assert(ixy==2);
	mu = 3;
	mv = 2;
	mw = 4;
    }

    const int n_rho = n_offset + _rho;
    const int nu    = n_offset + mu;
    const int nv    = n_offset + mv;
    const int nw    = n_offset + mw;
    const int n_N   = n_offset + _N;

    // Average states
    const double& gamma  = fiveMomentParams.get_gamma();
    const double& rho    = Q_ave.get(n_rho);
    const double  u1     = Q_ave.get(nu)/rho;
    const double  u2     = Q_ave.get(nv)/rho;
    const double  u3     = Q_ave.get(nw)/rho;
    const double& energy = Q_ave.get(n_N);
    const double  umag2  = (u1*u1 + u2*u2 + u3*u3);
    const double  press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
    const double  c      = sqrt(fabs(gamma*press/rho));
    const double  H      = (energy+press)/rho;  
    
    // Project onto right eigenvectors
    for (int k=1; k<=Q.getsize(2); k++)
    {
        Q.set(n_rho,k, W.get(n_rho,k) + W.get(nu,k) + W.get(n_N,k)  );

        Q.set(nu,k, (u1-c)*W.get(n_rho,k) + u1*W.get(nu,k) 
                  + (u1+c)*W.get(n_N,k) );

        Q.set(nv,k, u2*(W.get(n_rho,k) + W.get(nu,k) + W.get(n_N,k))
                  + W.get(nv,k) );
        
        Q.set(nw,k, u3*(W.get(n_rho,k) + W.get(nu,k) + W.get(n_N,k))
                  + W.get(nw,k) );

        Q.set(n_N, k, (H-u1*c)*W.get(n_rho,k) + umag2/2.0*W.get(nu,k)
                  + u2*W.get(nv,k) + u3*W.get(nw,k) 
                  + (H+u1*c)*W.get(n_N,k) );
    }
}

