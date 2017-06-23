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
// Shallow water equations.
//
// See also: $FINESS/lib/1d/ConstructIntegratedF.cpp.
void ConstructDiffTransformL( double dt, StateVars& Q, dTensorBC1& smax, dTensorBC2& F, dTensorBC2& Lstar)
{

    // Central difference routine (depends on spatial order!)
    // (One check to avoid if clause overhead)
    void (*CentralDifferences)( double dx, const dTensor2& f, dTensor2& fderivs) = GetCentralDifferences();

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
    const int MAX_DERIVS      = mpts_sten;         assert_ge( MAX_DERIVS, global_ini_params.get_space_order() );
    const int MAX_FLUX_DERIVS = mpts_sten-1;

    // Double check there are enough boundary points (otherwise there will be
    //                                                an unnoticed seg. fault)
    assert_ge( 1+mbc-mbc_small, half_mpts_sten );

    // Quadrature rules for numerically evaluating the integral of the flux
    dTensor1 w1d( MAX_DERIVS );
    dTensor1 x1d( MAX_DERIVS );

    void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d);
    void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);

    if( mpts_sten < 9 )
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

        // Save the function b_x:
        dTensor2 bxvals( 1, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                qvals.set( m, r, q.get( i+s, m ) );
                s++;
            }
        }
       if( maux >0)
       {
           int s = -half_mpts_sten+1;
           //printf("aux should be %e \n",aux.get(i+s,1));
           for( int r = 1; r <= mpts_sten; r++ )
           {
                bxvals.set( 1, r, aux.get( i+s, 1 ) );
                s++;
           }
        }
        else
        {
           int s = -half_mpts_sten+1;
           for( int r = 1; r <= mpts_sten; r++ )
           {
                bxvals.set( 1, r, 0.0 );
                s++;
           }
       
        }


        // Storage for local derivatives:
        dTensor2 qderivs ( meqn, MAX_DERIVS );
        
        // Storage for local for source:
        dTensor2 bxderivs ( 1, MAX_DERIVS );

        // Compute a FD approximation to the derivatives:
        CentralDifferences( dx, qvals, qderivs );

        // Compute a FD approximation to the derivatives of the bathymatry
        CentralDifferences( dx, bxvals, bxderivs );


        // Save all of the "zeroth" time derivatives.
        dTensor3 Q_mixed_derivs( meqn, MAX_DERIVS, MAX_DERIVS );
        dTensor2 BX_mixed_derivs(  MAX_DERIVS, MAX_DERIVS );BX_mixed_derivs.setall(0.0);
        dTensor3 G( 3, MAX_DERIVS, MAX_DERIVS );
        dTensor3 G1( 2, MAX_DERIVS, MAX_DERIVS );
 
        dTensor2 S( MAX_DERIVS, MAX_DERIVS );

        for( int me =1; me <= meqn; me++ )
        {
            for( int h=1; h <= MAX_DERIVS; h++ )
            { 
                Q_mixed_derivs.set( me, h, 1, qderivs.get(me, h)/factorial[h-1] ); 
            }
        }

        for( int h=1; h <= MAX_DERIVS; h++ )
        {
           BX_mixed_derivs.set( h, 1, bxderivs.get(1,h)/factorial[h-1] );
        }


        //Populate "zeroth" time derivatives of (phi*u)^2...
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                int k=0;
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 2, r+1, s+1)*Q_mixed_derivs.get(2, h+1-r, k+1-s );
                }
                G1.set( 1, h+1, 1, tmp );
        }

        //Populate "zeroth" time derivatives of h b_x
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                int k=0;
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, 1)*BX_mixed_derivs.get(h+1-r,k+1-s );
                }
                S.set( h+1, 1, tmp );
        }


        //Populate "zeroth" time derivatives of 1/u...
        G1.set(2,1,1,1.0/(Q_mixed_derivs.get(1,1,1)));
        for( int h=1; h < MAX_DERIVS; h++ )
        {
                int k=0;
                double tmp = 0.0;
                G1.set(2,h+1,1,0.0);
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*G1.get(2, h+1-r, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1));
                G1.set( 2, h+1, 1, tmp );
        }
       

        //compute (phi*u)^2/u and u^2
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                int k=0;
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                double tmp3 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += G1.get( 1, r+1, s+1)*G1.get(2, h+1-r, k+1-s );
                      tmp2 += Q_mixed_derivs.get( 1, r+1, s+1)*Q_mixed_derivs.get(1, h+1-r, k+1-s );
                }
                G.set(1,h+1,1,tmp1);
                G.set(2,h+1,1,tmp2);
                G.set(3,h+1,1,tmp3);       

         }
//assert_le( fabs( Q_mixed_derivs.get( 1, 1, 1 ) - q.get(i,1) ), 1e-13 );
//assert_le( fabs( Q_mixed_derivs.get( 2, 1, 1 ) - q.get(i,2) ), 1e-13 );

        // Recursive relationship goes here! //
        for( int k=0; k < MAX_FLUX_DERIVS;   k++ )      
        {

          //Populate time derivatives of (phi*u)^2...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 2, r+1, s+1)*Q_mixed_derivs.get(2, h+1-r, k+1-s );
                }
                G1.set( 1, h+1, k+1, tmp );
          }

          //Populate time derivatives of (h b_x)...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*BX_mixed_derivs.get(h+1-r, k+1-s );
                }
                S.set( h+1, k+1, tmp );
          }



          //Populate time derivatives of 1/u...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                G1.set(2,h+1,k+1,0.0);
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*G1.get(2, h+1-r, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1));
                G1.set( 2, h+1, k+1, tmp );
          }

          //compute (phi*u)^2/u and u^2
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                double tmp3 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += G1.get( 1, r+1, s+1)*G1.get(2, h+1-r, k+1-s );
                      tmp2 += Q_mixed_derivs.get( 1, r+1, s+1)*Q_mixed_derivs.get(1, h+1-r, k+1-s );
                }
                G.set(1,h+1,k+1,tmp1);
                G.set(2,h+1,k+1,tmp2);
                G.set(3,h+1,k+1,tmp3);

           }
          for( int h=0; h < MAX_FLUX_DERIVS-k; h++ )
          {

           //compute conserved quantity derivatives...
           Q_mixed_derivs.set(1,h+1,k+2,-((double)h+1.0)/(k+1)*Q_mixed_derivs.get(2,h+2,k+1) );
           Q_mixed_derivs.set(2,h+1,k+2,-((double)h+1.0)/(k+1)*(G.get(1,h+2,k+1)+G.get(2,h+2,k+1)+G.get(3,h+2,k+1)) - 1.0/((double)k+1.0)*S.get(h+1,k+1));
           }
         }


        // Construct the time-averaged flux.
        double u1_ta_flux = 0.;
        double u2_ta_flux = 0.;
        for( int mq=1; mq <= MAX_DERIVS; mq++ )
        {
            // Evaluate q at this quadrature point
            double u1 = Q_mixed_derivs.get(1,1,1);
            double u2 = Q_mixed_derivs.get(2,1,1);
            for( int k=1; k < MAX_FLUX_DERIVS; k++ )
            { 
                u1 += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 1, 1, k+1 ); 
                u2 += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 2, 1, k+1 ); 
            }

            // 0.5 * wgt * fluxfunc( q(t(xi)) ).
            u1_ta_flux += (0.5*w1d.get( mq )) * ( u2 );
            u2_ta_flux += (0.5*w1d.get( mq )) * ( 0.5*u1*u1+u2*u2/u1);
        }
        F.set( i, 1, u1_ta_flux );
        F.set( i, 2, u2_ta_flux );

        double source_ta = 0.;
        double source_ta_2 = 0.;

        for (int mq=1;mq<=MAX_DERIVS;mq++)
        {
            double u1 = S.get(1,1);
            double u2 = -S.get(1,1);
            for( int k=1; k < MAX_FLUX_DERIVS; k++ )
            {
                u1 +=  ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*S.get( 1, k+1 );
                u2 += -( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*S.get( 1, k+1 );
            }
            source_ta += (0.5*w1d.get( mq )) * u2;
        }
        Lstar.set(i,1,0.0);
        Lstar.set(i,2,source_ta);



    }

}

