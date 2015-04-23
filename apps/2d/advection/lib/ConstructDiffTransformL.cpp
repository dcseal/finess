#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "assert.h"
#include "IniParams.h"
#include "StateVars.h"
#include "CentralDifferences2D.h"

void ConstructDiffTransformL( double dt, StateVars& Q,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G)
{

    const dTensorBC3& q   = Q.const_ref_q  ();
    const dTensorBC3& aux = Q.const_ref_aux();

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double dy    = global_ini_params.get_dy();
    const double xlow  = global_ini_params.get_xlow();
    const double ylow  = global_ini_params.get_ylow();


    const int mpts_sten       = global_ini_params.get_space_order(); assert_eq( mpts_sten%2, 1 );

    const int mbc_small       = (mbc+1)/2;
    const int half_mpts_sten  = (mpts_sten+1)/2;          // assert_eq( half_mpts_sten, 3 );
    const int MAX_DERIVS      = mpts_sten;
    const int MAX_FLUX_DERIVS = mpts_sten-1;

    // Quadrature rules for numerically evaluating the integral of the flux
    const int NumQuad = 9;
    dTensor1 w1d( NumQuad );
    dTensor1 x1d( NumQuad );



const int ndim = 2;

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    for( int j = 1-mbc_small; j <= my+mbc_small; j++ )
    {

        // Physical location for this current value:
        dTensor2 xpts( 1, ndim );
        xpts.set( 1, 1, xlow + double(i)*dx - 0.5*dx );
        xpts.set( 1, 2, ylow + double(j)*dy - 0.5*dy );

        // Save the flux function:
        dTensor3 qvals ( meqn, mpts_sten, mpts_sten );
        //dTensor2 qvalsx ( meqn, mpts_sten );
        //dTensor2 qvalsy ( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s1 = -half_mpts_sten+1;
            for( int k = 1; k <= mpts_sten; k++ )
            {
               int s2 = -half_mpts_sten+1;
               for( int r = 1; r <= mpts_sten; r++ )
               {
                qvals.set( m, r,k, q.get( i+s2, j+s1, m ) );
//                printf("here %f %d %d %d \n",q.get( i+s2, j+s1, m ),i+s2,j+s1,m);
                s2++;
               }
               s1++;
            }
           
        }

        // Storage for local derivatives:
        dTensor3 qderivs ( meqn, MAX_DERIVS, MAX_DERIVS );

        // Compute a FD approximation to the derivatives:
        // Central difference routine (depends on spatial order!)
        void (*CentralDifferences2D)( double dx, double dy, const dTensor3& f, dTensor3& fderivs) = GetCentralDifferences2D();
        CentralDifferences2D( dx, dy, qvals, qderivs );
        // Save all of the "zeroth" time derivatives.
        dTensor4 Q_mixed_derivs( meqn, MAX_DERIVS, MAX_DERIVS, MAX_DERIVS );
        for( int me =1; me <= meqn; me++ )
        {
            for( int hx=1; hx <= MAX_DERIVS; hx++ )
            for( int hy=1; hy <= MAX_DERIVS; hy++ )
            {   double x1=xpts.get(1,1);double y1=xpts.get(1,1);double deriv1=1.0;
                Q_mixed_derivs.set( meqn, hx,hy, 1, qderivs.get(meqn, hx,hy)/(factorial[hx-1]*factorial[hy-1]) );
            }
        }


        // Recursive relationship goes here! //
        for( int me = 1; me <= meqn; me++ )
        {

            // TODO - THIS IS HARD CODED in two places: 
            //
            //    i) spatial (and hence temporal) order
            //   ii) velocity is set equal to one
            //
            for( int k=0; k < MAX_FLUX_DERIVS;   k++ )
            for( int hx=0; hx < MAX_FLUX_DERIVS-k; hx++ )
            for( int hy=0; hy < MAX_FLUX_DERIVS-k; hy++ )
            {
                double tmp = -(hx+1.0)/(k+1.0) * Q_mixed_derivs.get( me, hx+2, hy+1, k+1 ) -(hy+1.0)/(k+1.0) * Q_mixed_derivs.get( me, hx+1, hy+2, k+1 );
                Q_mixed_derivs.set( me, hx+1, hy+1, k+2, tmp );
            }
        }


        // Now, construct time-averaged flux.  Note that Q_mixed_derivs
        // already contains a factorial in its definition.  In order to
        // construct the "time-averaged" version, we need to only divide by a
        // single extra factor of (k+1).
        for( int m=1; m<=meqn; m++ )
        {
            double tmp = Q_mixed_derivs.get(m,1,1,1);
            for( int k=1; k < MAX_FLUX_DERIVS; k++ )
            {
                tmp += ( pow(dt,k) / (1.0+k) )*Q_mixed_derivs.get( m, 1, 1, k+1 );
//              tmp += (dt / 2.0) * Q_mixed_derivs.get( m, 1, k+1 );
            }
            F.set( i, j, m, tmp );
        }

        // Now, construct time-averaged flux.  Note that Q_mixed_derivs
        // already contains a factorial in its definition.  In order to
        // construct the "time-averaged" version, we need to only divide by a
        // single extra factor of (k+1).
        for( int m=1; m<=meqn; m++ )
        {
            double tmp = Q_mixed_derivs.get(m,1,1,1);
            for( int k=1; k < MAX_FLUX_DERIVS; k++ )
            {
                tmp += ( pow(dt,k) / (1.0+k) )*Q_mixed_derivs.get( m, 1, 1, k+1 );
//              tmp += (dt / 2.0) * Q_mixed_derivs.get( m, 1, k+1 );
            }
            G.set( i, j, m, tmp );
        }


        }


}


