#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"
#include "StateVars.h"
#include "assert.h"

// Time expanded state variable, q using discrete transform.
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

// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );
const int MAX_DERIVS     = 5;
const int MAX_FLUX_DERIVS = 4;

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
        void CentralDifferences( double dx, const dTensor2& f, dTensor2& fderivs);
        CentralDifferences( dx, qvals, qderivs );

        // Save all of the "zeroth" time derivatives.
        dTensor3 Q_mixed_derivs( meqn, MAX_DERIVS, MAX_DERIVS );
        for( int me =1; me <= meqn; me++ )
        {
            for( int h=1; h <= MAX_DERIVS; h++ )
            { 
                Q_mixed_derivs.set( me, h, 1, qderivs.get(me, h)/factorial[h-1] ); 
            }
        }

//assert_le( fabs( Q_mixed_derivs.get( 1, 1, 1 ) - q.get(i,1) ), 1e-13 );
//assert_le( fabs( Q_mixed_derivs.get( 2, 1, 1 ) - q.get(i,2) ), 1e-13 );

        // Easy index for the two-variable system.
        int b4[2]={2,1};

        // Recursive relationship goes here! //
        for( int k=0; k < MAX_FLUX_DERIVS;   k++ )      
        {
           for( int me = 1; me <= meqn; me++ )
           {

               for( int h=0; h < MAX_FLUX_DERIVS-k; h++ )
               {
                  double tmp = -(h+1.0)/(k+1.0) * Q_mixed_derivs.get( b4[me-1], h+2, k+1 );
                  Q_mixed_derivs.set( me, h+1, k+2, tmp );
               }
           }
         }

        // Now, construct time-averaged flux.  Note that Q_mixed_derivs
        // already contains a factorial in its definition.  In order to
        // construct the "time-averaged" version, we need to only divide by a
        // single extra factor of (k+1).
        for( int m=1; m<=meqn; m++ )
        {
            double tmp = Q_mixed_derivs.get(m,1,1);
            for( int k=1; k < MAX_FLUX_DERIVS; k++ )
            {
                tmp += ( pow(dt,k) / (1.0+k) )*Q_mixed_derivs.get( m, 1, k+1 );
//              tmp += (dt / 2.0) * Q_mixed_derivs.get( m, 1, k+1 );
            }

            // Because we are taking Taylor expansions of the *fluxes*, and
            // not the conserved variables, we need to "swap" these two one
            // more time.  That is, f^n = ( q2, q1 ), not (q1, q2 ).
//          F.set( i, m, tmp );
            F.set( i, b4[m-1], tmp );

        }

    }

}

