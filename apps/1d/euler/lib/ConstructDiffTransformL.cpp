#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"
#include "StateVars.h"
#include "assert.h"
#include "CentralDifferences.h"

// Time expanded state variable, q using discrete transform.
// See: Multi-moment ADER-Taylor methods for systems of conservation laws
// with source terms in one dimension 
//
// Euler equations.
//
// See also: $FINESS/lib/1d/ConstructIntegratedF.cpp.
void ConstructDiffTransformL( double dt, StateVars& Q, dTensorBC1& smax, dTensorBC2& F)
{

    dTensorBC2& q   = Q.ref_q();
    dTensorBC2& aux = Q.ref_aux();

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    const int mpts_sten       = global_ini_params.get_space_order(); assert_eq( mpts_sten%2, 1 );

    const int mbc_small       = (mbc+1)/2;
    const int half_mpts_sten  = (mpts_sten+1)/2;          // assert_eq( half_mpts_sten, 3 );
    const int MAX_DERIVS      = mpts_sten;
    const int MAX_FLUX_DERIVS = mpts_sten-1;

    // Quadrature rules for numerically evaluating the integral of the flux
    const int NumQuad = 9;
    dTensor1 w1d( NumQuad );
    dTensor1 x1d( NumQuad );

    void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d);
    void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);

    //if( mpts_sten < 9 )
    if( NumQuad < 9 )
        setGaussLobattoPoints1d( w1d, x1d);
    else
        setGaussPoints1d( w1d, x1d );

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    {

        // Physical location for this current value:
        dTensor1 xpts( 1 );
        xpts.set( 1, xlow + double(i)*dx - 0.5*dx );

        // Save the flux function:
        dTensor2 qvals( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                qvals.set( m, r, q.get( i+s, m ) );
                s++;
            }
        }

        // Storage for local derivatives:
        dTensor2 qderivs ( meqn, MAX_DERIVS );

        // Compute a FD approximation to the derivatives:
        // Central difference routine (depends on spatial order!)
        void (*CentralDifferences)( double dx, const dTensor2& f, dTensor2& fderivs) = GetCentralDifferences();
        CentralDifferences( dx, qvals, qderivs );

        // Save all of the "zeroth" time derivatives.
        dTensor3 Q_mixed_derivs( meqn, MAX_DERIVS, MAX_DERIVS );
        dTensor3 Gs( 1, MAX_DERIVS, MAX_DERIVS );
        dTensor3 R( 1, MAX_DERIVS, MAX_DERIVS );
        dTensor3 Gr( 2, MAX_DERIVS, MAX_DERIVS );
        dTensor3 Ge( 1, MAX_DERIVS, MAX_DERIVS );
        dTensor3 P( 1, MAX_DERIVS, MAX_DERIVS );

        //Initialize all arrays to zero(maybe done already...i should check the constructor
        for( int h=1; h <= MAX_DERIVS; h++ )
        for( int k=1; k <= MAX_DERIVS; k++ )
        {
            for( int me =1; me <= meqn; me++ )
            {
                Q_mixed_derivs.set( me, h, k, 0.0 );
            }
            Gs.set(1,h,k,0.0);
            R.set(1,h,k,0.0);
            Gr.set(1,h,k,0.0);
            Gr.set(2,h,k,0.0);
            Ge.set(1,h,k,0.0);
            P.set(1,h,k,0.0);
        }
        for( int me =1; me <= meqn; me++ )
        {
            for( int h=1; h <= MAX_DERIVS; h++ )
            { 
                Q_mixed_derivs.set( me, h, 1, qderivs.get(me, h)/factorial[h-1] ); 
            }
        }

//assert_le( fabs( Q_mixed_derivs.get( 1, 1, 1 ) - q.get(i,1) ), 1e-13 );
//assert_le( fabs( Q_mixed_derivs.get( 2, 1, 1 ) - q.get(i,2) ), 1e-13 );

        //Need to deal with k=0 specially because of the 1/rho term...
        for(int k=0;k<1;k++)
        {


          //Populate time derivatives of (rho u)^2...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 2, r+1, s+1)*Q_mixed_derivs.get(2, h+1-r, k+1-s );
                }
                Gs.set( 1, h+1, k+1, tmp );
          }

          R.set( 1, 1, 1, 1.0/(Q_mixed_derivs.get(1,1,1)) );
          //Populate time derivatives of 1/rho...
          for( int h=1; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                R.set(1,h+1,k+1,0.0);
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1));
                R.set( 1, h+1, k+1, tmp );
          }
          //populate the time derivatives of (rho u)^2/rho and u
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Gs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                      tmp1 += Q_mixed_derivs.get( 2, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                Gr.set( 1, h+1, k+1, tmp );
                Gr.set( 2, h+1, k+1, tmp1 );
          }  
          //Populate time derivatives of pressure
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                tmp = 0.4*(Q_mixed_derivs.get( 5, h+1, k+1)-0.5*Gr.get(1,h+1,k+1));
                P.set( 1, h+1, k+1, tmp );
          }
       
          //populate time derivatives of energy equation flux...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += (Q_mixed_derivs.get( 5, r+1, s+1)+P.get( 1, r+1, s+1))*Gr.get(2, h+1-r, k+1-s );
                }
                Ge.set(1,h+1,k+1,tmp1);

           }
          for( int h=0; h < MAX_FLUX_DERIVS-k; h++ )
          {

           //compute conserved quantity derivatives...
           Q_mixed_derivs.set(1,h+1,k+2,-((double)h+1.0)/(k+1)*Q_mixed_derivs.get(2,h+2,k+1));
           Q_mixed_derivs.set(2,h+1,k+2,-((double)h+1.0)/(k+1)*(Gr.get(1,h+2,k+1)+P.get(1,h+2,k+1)));
           Q_mixed_derivs.set(3,h+1,k+2,0.0);
           Q_mixed_derivs.set(4,h+1,k+2,0.0);
           Q_mixed_derivs.set(5,h+1,k+2,-((double)h+1.0)/(k+1)*Ge.get(1,h+2,k+1));
           }

       }


        // Recursive relationship goes here! Compute higher order derivatives //
        for( int k=1; k < MAX_FLUX_DERIVS;   k++ )      
        {
        //Populate time derivatives of (rho u)^2...
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 2, r+1, s+1)*Q_mixed_derivs.get(2, h+1-r, k+1-s );
                }
                Gs.set( 1, h+1, k+1, tmp );
        }
        //Populate time derivatives of 1/rho...
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp = 0.0;
                R.set(1,h+1,k+1,0.0);
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1));
                R.set( 1, h+1, k+1, tmp );
        }

        //populate time derivatives of (rho u)^2/rho and u
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp = 0.0;
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Gs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                      tmp1 += Q_mixed_derivs.get( 2, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                Gr.set( 1, h+1, k+1, tmp );
                Gr.set( 2, h+1, k+1, tmp1 );
        }
        //populate time derivatives of P
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp = 0.0;
                tmp = 0.4*(Q_mixed_derivs.get( 5, h+1, k+1)-0.5*Gr.get(1,h+1,k+1));
                P.set( 1, h+1, k+1, tmp );
        }
        //populate time derivatives of energy equation flux
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += (Q_mixed_derivs.get( 5, r+1, s+1)+P.get( 1, r+1, s+1))*Gr.get(2, h+1-r, k+1-s );
                }
                Ge.set(1,h+1,k+1,tmp1);

         }

          for( int h=0; h < MAX_FLUX_DERIVS-k; h++ )
          {

           //compute conserved quantity derivatives...
           Q_mixed_derivs.set(1,h+1,k+2,-((double)h+1.0)/(k+1)*Q_mixed_derivs.get(2,h+2,k+1));
           Q_mixed_derivs.set(2,h+1,k+2,-((double)h+1.0)/(k+1)*(Gr.get(1,h+2,k+1)+P.get(1,h+2,k+1)));
           Q_mixed_derivs.set(3,h+1,k+2,0.0);
           Q_mixed_derivs.set(4,h+1,k+2,0.0);
           Q_mixed_derivs.set(5,h+1,k+2,-((double)h+1.0)/(k+1)*Ge.get(1,h+2,k+1));
           }
         }
        //Compute the highest time-order terms of the auxillary quantities because we will need these for the Euler fluxes 
        for(int k=MAX_FLUX_DERIVS-1;k<MAX_DERIVS;k++)
        {


          //Populate time derivatives of (rho u)^2...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 2, r+1, s+1)*Q_mixed_derivs.get(2, h+1-r, k+1-s );
                }
                Gs.set( 1, h+1, k+1, tmp );
          }
          //Populate time derivatives of 1/rho...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                R.set(1,h+1,k+1,0.0);
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1));
                R.set( 1, h+1, k+1, tmp );
          }
          //populate time derivatives of (rho u)^2/rho and u
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Gs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                      tmp1 += Q_mixed_derivs.get( 2, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                Gr.set( 1, h+1, k+1, tmp );
                Gr.set( 2, h+1, k+1, tmp1 );
          }  
          //populate time derivatives of pressure
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                tmp = 0.4*(Q_mixed_derivs.get( 5, h+1, k+1)-0.5*Gr.get(1,h+1,k+1));
                P.set( 1, h+1, k+1, tmp );
          }
       
          //populate time derivatives of energy equation flux function
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += (Q_mixed_derivs.get( 5, r+1, s+1)+P.get( 1, r+1, s+1))*Gr.get(2, h+1-r, k+1-s );
                }
                Ge.set(1,h+1,k+1,tmp1);

           }

       }


        // Construct the time-averaged flux.
        double rho_ta_flux = 0.;
        double u1_ta_flux = 0.;
        double E_ta_flux = 0.;
        for( int mq=1; mq <= NumQuad; mq++ )
        {
            // Evaluate q at this quadrature point
            double rho= Q_mixed_derivs.get(1,1,1);
            double ru1 = Q_mixed_derivs.get(2,1,1);
            double E = Q_mixed_derivs.get(5,1,1);
            double P1 = P.get(1,1,1);
            double Ge1 = Ge.get(1,1,1);
            for( int k=1; k < MAX_DERIVS; k++ )
            { 
                rho += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 1, 1, k+1 ); 
                ru1 += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 2, 1, k+1 ); 
                E += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 5, 1, k+1 );
                //the pressure function 
                P1 += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*P.get( 1, 1, k+1 ); 
                //the function that gives the flux for the energy equation
                Ge1 += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Ge.get( 1, 1, k+1 ); 
            }
            rho_ta_flux += (0.5*w1d.get( mq )) * ( ru1 );
            u1_ta_flux += (0.5*w1d.get( mq )) * ( ru1*ru1/rho+P1);
            E_ta_flux += (0.5*w1d.get( mq )) * Ge1;
        }
        F.set( i, 1, rho_ta_flux );
        F.set( i, 2, u1_ta_flux );
        F.set( i, 3, 0.0 );
        F.set( i, 4, 0.0 );
        F.set( i, 5, E_ta_flux );

    }

}

