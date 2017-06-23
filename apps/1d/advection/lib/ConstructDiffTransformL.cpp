#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"
#include "StateVars.h"
#include "assert.h"
#include "CentralDifferences.h"

using namespace std;

// Time expanded state variable, q using discrete transform.
//
// Advection equation, q_t + ( a q )_x = 0.
//
// This routine only solves the constant coefficient case, and it assumes that
// the advection velocity is stored in aux(1,1,1).
//
// See also: $FINESS/lib/1d/ConstructIntegratedF.cpp.
//void ConstructDiffTransformL( double dt, StateVars& Q, dTensorBC1& smax, dTensorBC2& F)
void ConstructDiffTransformL( double dt, StateVars& Q, dTensorBC1& smax, dTensorBC2& F, dTensorBC2& Lstar)
{

    // Central difference routine (depends on spatial order!)
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
    const int MAX_DERIVS      = mpts_sten;
    const int MAX_FLUX_DERIVS = mpts_sten-1;

    // Quadrature rules for numerically evaluating the integral of the flux
    dTensor1 w1d( MAX_DERIVS );
    dTensor1 x1d( MAX_DERIVS );

    void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d);
    void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);

    if( mpts_sten < 9 )
        setGaussLobattoPoints1d( w1d, x1d);
    else
        setGaussPoints1d( w1d, x1d );


    dTensorBC2 flag(mx,meqn,mbc);
    void FlagIndicator( StateVars& Q, dTensorBC2& flag );
    FlagIndicator( Q, flag );

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
        CentralDifferences( dx, qvals, qderivs );

        // Save all of the "zeroth" time derivatives.
        dTensor3 Q_mixed_derivs( meqn, MAX_DERIVS, MAX_DERIVS );
        for( int me =1; me <= meqn; me++ )
        {
            for( int h=1; h <= MAX_DERIVS; h++ )
            { 
                Q_mixed_derivs.set( meqn, h, 1, qderivs.get(meqn, h)/factorial[h-1] ); 
            }
        }

// assert_le( fabs( Q_mixed_derivs.get( 1, 1, 1 ) - q.get(i,1) ), 1e-13 );

        // Recursive relationship goes here! //
        for( int me = 1; me <= meqn; me++ )
        {

            // TODO - THIS IS HARD CODED in two places: 
            //
            //    i) spatial (and hence temporal) order
            //   ii) velocity is set equal to one
            //
            for( int k=0; k < MAX_FLUX_DERIVS;   k++ )      
            for( int h=0; h < MAX_FLUX_DERIVS-k; h++ )
            {
                double tmp = -(h+1.0)/(k+1.0) * Q_mixed_derivs.get( me, h+2, k+1 );
                Q_mixed_derivs.set( me, h+1, k+2, tmp );
            }
        }

        // Now, construct time-averaged flux.  Note that Q_mixed_derivs
        // already contains a factorial in its definition.  In order to
        // construct the "time-averaged" version, we need to only divide by a
        // single extra factor of (k+1).
        for( int m=1; m<=meqn; m++ )
        {
            const int nterms = 4*flag.get(i,m) + MAX_FLUX_DERIVS*(1-flag.get(i,m));
            double tmp = Q_mixed_derivs.get(m,1,1);
//          for( int k=1; k < MAX_FLUX_DERIVS; k++ )
            for( int k=1; k < nterms; k++ )
            {
                tmp += ( pow(dt,k) / (1.0+k) )*Q_mixed_derivs.get( m, 1, k+1 );
            }
            F.set( i, m, tmp );
        }

        // Construct the time-averaged flux.
//      double ta_flux = 0.;
//      for( int m=1; m<=meqn; m++ )
//      {
//          for( int mq=1; mq <= MAX_DERIVS; mq++ )
//          {
//              // Evaluate q at this quadrature point
//              double q = Q_mixed_derivs.get(m,1,1);
//              for( int k=1; k < MAX_FLUX_DERIVS; k++ )
//              { q += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( m, 1, k+1 ); }

//              // 0.5 * wgt * fluxfunc( q(t(xi)) ).
//              ta_flux += (0.5*w1d.get( mq )) * ( q );
//          }
//          F.set( i, m, ta_flux );
//      }

    }

}

void FlagIndicator( StateVars& Q, dTensorBC2& flag )
{

    void SetBndValues( StateVars& Q );
    SetBndValues( Q );
    dTensorBC2& q   = Q.ref_q();
    dTensorBC2& aux = Q.ref_aux();

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    // Multiresolution analysis parameter epsilon
    const double eps_mr  = 0.1;

    dTensorBC2 flag_no_buff( mx, meqn, mbc );
    flag_no_buff.setall(0.);

    // --------------------------------------------------------------------
    // Part I: Fill in values for the flag
    // --------------------------------------------------------------------
    int num_flagged_elems = 0;
    for( int i = 1-mbc+1; i <= mx+mbc-1; i++ )
    {

        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i+1,m) + q.get(i-1,m) );
            Qavg.set(m, tmp );

            if( abs( tmp - q.get(i,m) ) > eps_mr * dx )
            { 
                flag_no_buff.set(i, m, 1); 
                num_flagged_elems += 1;     
//              cout << "flagged elem i  = " << i << endl;
            }

        }

    }
    cout << "num_flagged_elems = " << num_flagged_elems << endl;

    // --------------------------------------------------------------------
    // Part II: Create buffer zone around each flagged element
    // --------------------------------------------------------------------
    flag.setall(0.);
    for( int i = 1; i <= mx; i++ )
    {

        for( int m=1; m <= meqn; m++ )
        {
            if( flag_no_buff.get(i,m) > 0 )
            { 
                for( int j=-mbc; j <= mbc; j++ )
                {
                    flag.set(i+j, m, 1); 
                }
            }

        }

    }

}
