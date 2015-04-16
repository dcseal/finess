#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"
#include "StateVars.h"
#include "assert.h"
#include "CentralDifferences.h"

// Time expanded state variable, q using discrete transform.
//
// Burgers equation, q_t + (0.5 q^2)_x = 0.
//
// See also: $FINESS/lib/1d/ConstructIntegratedF.cpp.
void ConstructDiffTransformL( double dt, StateVars& Q, dTensorBC1& smax, dTensorBC2& F)
{

    // Central difference routine (depends on spatial order!)
    void (*CentralDifferences)( double dx, const dTensor2& f, dTensor2& fderivs) = GetCentralDifferences();

    dTensorBC2& q   = Q.ref_q();
    dTensorBC2& aux = Q.ref_aux();

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);  assert_eq( meqn, 1 );  const int me = 1;
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

/*
// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );
const int MAX_DERIVS     = 5;
const int MAX_FLUX_DERIVS = 4;
*/

    // Quadrature rules for numerically evaluating the integral of the flux
    dTensor1 w1d( MAX_DERIVS );
    dTensor1 x1d( MAX_DERIVS );
    void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d);
    setGaussLobattoPoints1d( w1d, x1d);

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    {

        // Physical location for this current value:
        dTensor1 xpts( 1 );
        xpts.set( 1, xlow + double(i)*dx - 0.5*dx );

        // Save the flux function:
        dTensor2 qvals( me, mpts_sten );

        int s = -half_mpts_sten+1;
        for( int r = 1; r <= mpts_sten; r++ )
        {
            qvals.set( me, r, q.get( i+s, me ) );
            s++;
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

        // Recursive relationship goes here! //
        for( int k=0; k < MAX_FLUX_DERIVS;   k++ )      
        for( int h=0; h < MAX_FLUX_DERIVS-k; h++ )
        {

            // Compute discrete transform of q*q:
            double tmp = 0.0;
            for( int r=0; r <= h+1; r++ )
            for( int s=0; s <= k; s++ )
            {
                tmp += Q_mixed_derivs.get( me, r+1, s+1)*Q_mixed_derivs.get(me, h+2-r, k+1-s );
            }
            tmp *= -0.5;
            Q_mixed_derivs.set( me, h+1, k+2, tmp );
        }

        // Construct the time-averaged flux.
        double ta_flux = 0.;
        for( int mq=1; mq <= MAX_DERIVS; mq++ )
        {
            // Evaluate q at this quadrature point
            double q = Q_mixed_derivs.get(me,1,1);
            for( int k=1; k < MAX_FLUX_DERIVS; k++ )
            { q += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( me, 1, k+1 ); }

            // 0.5 * wgt * fluxfunc( q(t(xi)) ).
            ta_flux += (0.5*w1d.get( mq )) * ( 0.5*q*q );
        }
        F.set( i, me, ta_flux );

    }

}

