#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "assert.h"
#include "IniParams.h"
#include "StateVars.h"
#include "CentralDifferences2D.h"

// Time expanded state variable of the flux using discrete differential 
// transforms
//
// 2D Euler equations.
//
// See also: $FINESS/lib/2d/ConstructIntegratedF.cpp.
void ConstructDiffTransformL( double dt, StateVars& Q, dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G)
{

    // Number of dimensions
    const int ndim = 2;

    // Central difference routine (depends on spatial order!)
    void (*CentralDifferences2D)( double dx, double dy, const dTensor3& f, dTensor3& fderivs) = GetCentralDifferences2D();

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
    double const gamma = global_ini_params.get_gamma();

    const int mpts_sten       = global_ini_params.get_space_order(); assert_eq( mpts_sten%2, 1 );

    const int mbc_small       = (mbc+1)/2;
    const int half_mpts_sten  = (mpts_sten+1)/2;          // assert_eq( half_mpts_sten, 3 );
    const int MAX_DERIVS      = mpts_sten;
    const int MAX_FLUX_DERIVS = mpts_sten-1;

    // Double check there are enough boundary points (otherwise there will be
    //                                                an unnoticed seg. fault)
    assert_ge( 1+mbc-mbc_small, half_mpts_sten );

    // Quadrature rules for numerically evaluating the integral of the flux
    const int NumQuad = 9;
    dTensor1 w1d( NumQuad );
    dTensor1 x1d( NumQuad );

    // Extra variables for nondimensional version of PDE.  This helps to
    // reduce the total number of divisions that happen inside the main loop
    const double dt_dx = dt/dx;
    const double dt_dy = dt/dy;

    void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d);
    void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);
    if( NumQuad < 9 )
        setGaussLobattoPoints1d( w1d, x1d);
    else
        setGaussPoints1d( w1d, x1d );

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
        qvals.setall(0.0);

        for( int m=1; m <= meqn; m++ )
        {
            int s1 = -half_mpts_sten+1;
            for( int k = 1; k <= mpts_sten; k++ )
            {
               int s2 = -half_mpts_sten+1;
               for( int r = 1; r <= mpts_sten; r++ )
               {
                   qvals.set( m, r,k, q.get( i+s2, j+s1, m ) );
                   // printf("here %f %d %d %d \n",q.get( i+s2, j+s1, m ),i+s2,j+s1,m);
                   s2++;
               }
               s1++;
            }
           
        }

        // Storage for local derivatives:
        dTensor3 qderivs ( meqn, MAX_DERIVS, MAX_DERIVS );
        qderivs.setall(0.0);

        // Compute a FD approximation to the derivatives:
        CentralDifferences2D( dx, dy, qvals, qderivs );

        // These need to be allocated here in order to be thread safe
        dTensor4 Gs  ( 3, MAX_DERIVS, MAX_DERIVS, MAX_DERIVS );
        dTensor4 R   ( 1, MAX_DERIVS, MAX_DERIVS, MAX_DERIVS );
        dTensor4 G2r ( 3, MAX_DERIVS, MAX_DERIVS, MAX_DERIVS );
        dTensor4 G1r ( 2, MAX_DERIVS, MAX_DERIVS, MAX_DERIVS );
        dTensor4 Ge  ( 2, MAX_DERIVS, MAX_DERIVS, MAX_DERIVS );
        dTensor4 P   ( 1, MAX_DERIVS, MAX_DERIVS, MAX_DERIVS );

        // Initialize all arrays to zero (maybe done already...i should check the constructor)
        Gs.setall ( 0. );
        R.setall  ( 0. );
        G2r.setall( 0. );
        G1r.setall( 0. );
        Ge.setall ( 0. );
        P.setall  ( 0. );


        // Save all of the "zeroth" time derivatives.
        dTensor4 Q_mixed_derivs( meqn, MAX_DERIVS, MAX_DERIVS, MAX_DERIVS );
        Q_mixed_derivs.setall(0.);
        for( int me =1; me <= meqn; me++ )
        {
            for( int hx=1; hx <= MAX_DERIVS; hx++ )
            for( int hy=1; hy <= MAX_DERIVS; hy++ )
            {   
                double x1     = xpts.get(1,1);
                double y1     = xpts.get(1,1);
                double deriv1 = 1.0;

                //Q_mixed_derivs.set( me, hx,hy, 1,qderivs.get(me, hx,hy)/(factorial[hx-1]*factorial[hy-1]) );
                //No need to do anything here because we have pulled division by the factorial into the central differences code.
                Q_mixed_derivs.set( me, hx, hy, 1, qderivs.get(me,hx,hy) );

            }
        }

        //Need to deal with k=0 specially because of the 1/rho term...
        for(int k=0;k<1;k++)
        {

          //Populate time derivatives of (rho u)^2...
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                double tmp3 = 0.0;
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*Q_mixed_derivs.get(2, hx+1-rx,hy+1-ry, k+1-s );
                      tmp2 += Q_mixed_derivs.get( 3, rx+1,ry+1, s+1)*Q_mixed_derivs.get(3, hx+1-rx,hy+1-ry, k+1-s );
                      tmp3 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*Q_mixed_derivs.get(3, hx+1-rx,hy+1-ry, k+1-s );
                }
                Gs.set( 1, hx+1,hy+1, k+1, tmp1 );
                Gs.set( 2, hx+1,hy+1, k+1, tmp2 );
                Gs.set( 3, hx+1,hy+1, k+1, tmp3 );
          }

          R.set( 1, 1,1, 1, 1.0/(Q_mixed_derivs.get(1,1,1,1)) );

          // Populate time derivatives of 1/rho...
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                if(hx+hy+k>0)
                {
                    double tmp = 0.0;
                    R.set(1,hx+1,hy+1,k+1,0.0);
                    for( int rx=0; rx <= hx; rx++ )
                    for( int ry=0; ry <= hy; ry++ )
                    for( int s=0; s <= k; s++ )
                    {
                        tmp += Q_mixed_derivs.get( 1, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                    }
                    tmp *=-1.0/(Q_mixed_derivs.get(1,1,1,1));
                    R.set( 1, hx+1,hy+1, k+1, tmp );
                }
          }

          // Populate the time derivatives of (rho u)^2/rho and u
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp1  = 0.0;
                double tmp2  = 0.0;
                double tmp3  = 0.0;
                double tmp11 = 0.0;
                double tmp12 = 0.0;
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1  += Gs.get( 1, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp2  += Gs.get( 2, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp3  += Gs.get( 3, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp11 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                      tmp12 += Q_mixed_derivs.get( 3, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                }
                G2r.set( 1, hx+1,hy+1, k+1, tmp1 );
                G2r.set( 2, hx+1,hy+1, k+1, tmp2 );
                G2r.set( 3, hx+1,hy+1, k+1, tmp3 );
                G1r.set( 1, hx+1,hy+1, k+1, tmp11 );
                G1r.set( 2, hx+1,hy+1, k+1, tmp12 );
          }

          //Populate time derivatives of pressure
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp = 0.0;
                tmp = (gamma-1.0)*(Q_mixed_derivs.get( 5, hx+1,hy+1, k+1)-0.5*G2r.get(1,hx+1,hy+1,k+1)-0.5*G2r.get(2,hx+1,hy+1,k+1));
                P.set( 1, hx+1,hy+1, k+1, tmp );
          }
       
          // Populate time derivatives of energy equation flux...
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += (Q_mixed_derivs.get( 5, rx+1,ry+1, s+1)+P.get( 1, rx+1,ry+1, s+1))*G1r.get(1, hx+1-rx,hy+1-ry, k+1-s );
                      tmp2 += (Q_mixed_derivs.get( 5, rx+1,ry+1, s+1)+P.get( 1, rx+1,ry+1, s+1))*G1r.get(2, hx+1-rx,hy+1-ry, k+1-s );
                }
                Ge.set(1,hx+1,hy+1,k+1,tmp1);
                Ge.set(2,hx+1,hy+1,k+1,tmp2);

           }

           for( int hy=0; hy < MAX_FLUX_DERIVS-k; hy++ )
           for( int hx=0; hx < MAX_FLUX_DERIVS-k; hx++ )
           {

               //compute conserved quantity derivatives...
               Q_mixed_derivs.set(1,hx+1,hy+1,k+2,-dt_dx*((double)hx+1.0)/(k+1)*Q_mixed_derivs.get(2,hx+2,hy+1,k+1)              -dt_dy*((double)hy+1.0)/(k+1)*Q_mixed_derivs.get(3,hx+1,hy+2,k+1));
               Q_mixed_derivs.set(2,hx+1,hy+1,k+2,-dt_dx*((double)hx+1.0)/(k+1)*(G2r.get(1,hx+2,hy+1,k+1)+P.get(1,hx+2,hy+1,k+1))-dt_dy*((double)hy+1.0)/(k+1)*G2r.get(3,hx+1,hy+2,k+1));
               Q_mixed_derivs.set(3,hx+1,hy+1,k+2,-dt_dy*((double)hy+1.0)/(k+1)*(G2r.get(2,hx+1,hy+2,k+1)+P.get(1,hx+1,hy+2,k+1))-dt_dx*((double)hx+1.0)/(k+1)*G2r.get(3,hx+2,hy+1,k+1));
               Q_mixed_derivs.set(4,hx+1,hy+1,k+2,0.0);
               Q_mixed_derivs.set(5,hx+1,hy+1,k+2,-dt_dx*((double)hx+1.0)/(k+1)*Ge.get(1,hx+2,hy+1,k+1)-dt_dy*((double)hy+1.0)/(k+1)*Ge.get(2,hx+1,hy+2,k+1));
           }

       }

       for(int k=1;k<MAX_FLUX_DERIVS;k++)
       {

           //Populate time derivatives of (rho u)^2...
           for( int hx=0; hx < MAX_DERIVS; hx++ )
           for( int hy=0; hy < MAX_DERIVS; hy++ )
           {
               double tmp1 = 0.0;
               double tmp2 = 0.0;
               double tmp3 = 0.0;
               for( int rx=0; rx <= hx; rx++ )
               for( int ry=0; ry <= hy; ry++ )
               for( int s=0; s <= k; s++ )
               {
                   tmp1 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*Q_mixed_derivs.get(2, hx+1-rx,hy+1-ry, k+1-s );
                   tmp2 += Q_mixed_derivs.get( 3, rx+1,ry+1, s+1)*Q_mixed_derivs.get(3, hx+1-rx,hy+1-ry, k+1-s );
                   tmp3 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*Q_mixed_derivs.get(3, hx+1-rx,hy+1-ry, k+1-s );
               }
               Gs.set( 1, hx+1,hy+1, k+1, tmp1 );
               Gs.set( 2, hx+1,hy+1, k+1, tmp2 );
               Gs.set( 3, hx+1,hy+1, k+1, tmp3 );
          }

          //Populate time derivatives of 1/rho...
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                //if(hx+hy>0)
                {

                double tmp = 0.0;
                R.set(1,hx+1,hy+1,k+1,0.0);
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1,1));
                R.set( 1, hx+1,hy+1, k+1, tmp );
                }
          }

          //populate the time derivatives of (rho u)^2/rho and u
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                double tmp3 = 0.0;
                double tmp11 = 0.0;
                double tmp12 = 0.0;
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += Gs.get( 1, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp2 += Gs.get( 2, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp3 += Gs.get( 3, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp11 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                      tmp12 += Q_mixed_derivs.get( 3, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                }
                G2r.set( 1, hx+1,hy+1, k+1, tmp1 );
                G2r.set( 2, hx+1,hy+1, k+1, tmp2 );
                G2r.set( 3, hx+1,hy+1, k+1, tmp3 );
                G1r.set( 1, hx+1,hy+1, k+1, tmp11 );
                G1r.set( 2, hx+1,hy+1, k+1, tmp12 );
          }

          //Populate time derivatives of pressure
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp = 0.0;
                tmp = (gamma-1.0)*(Q_mixed_derivs.get( 5, hx+1,hy+1, k+1)-0.5*G2r.get(1,hx+1,hy+1,k+1)-0.5*G2r.get(2,hx+1,hy+1,k+1));
                P.set( 1, hx+1,hy+1, k+1, tmp );
          }
       
          //populate time derivatives of energy equation flux...
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += (Q_mixed_derivs.get( 5, rx+1,ry+1, s+1)+P.get( 1, rx+1,ry+1, s+1))*G1r.get(1, hx+1-rx,hy+1-ry, k+1-s );
                      tmp2 += (Q_mixed_derivs.get( 5, rx+1,ry+1, s+1)+P.get( 1, rx+1,ry+1, s+1))*G1r.get(2, hx+1-rx,hy+1-ry, k+1-s );
                }
                Ge.set(1,hx+1,hy+1,k+1,tmp1);
                Ge.set(2,hx+1,hy+1,k+1,tmp2);

           }
          for( int hy=0; hy < MAX_FLUX_DERIVS-k; hy++ )
          for( int hx=0; hx < MAX_FLUX_DERIVS-k; hx++ )
          {

               //compute conserved quantity derivatives...
               Q_mixed_derivs.set(1,hx+1,hy+1,k+2,-dt_dx*((double)hx+1.0)/(k+1)*Q_mixed_derivs.get(2,hx+2,hy+1,k+1)              -dt_dy*((double)hy+1.0)/(k+1)*Q_mixed_derivs.get(3,hx+1,hy+2,k+1));
               Q_mixed_derivs.set(2,hx+1,hy+1,k+2,-dt_dx*((double)hx+1.0)/(k+1)*(G2r.get(1,hx+2,hy+1,k+1)+P.get(1,hx+2,hy+1,k+1))-dt_dy*((double)hy+1.0)/(k+1)*G2r.get(3,hx+1,hy+2,k+1));
               Q_mixed_derivs.set(3,hx+1,hy+1,k+2,-dt_dy*((double)hy+1.0)/(k+1)*(G2r.get(2,hx+1,hy+2,k+1)+P.get(1,hx+1,hy+2,k+1))-dt_dx*((double)hx+1.0)/(k+1)*G2r.get(3,hx+2,hy+1,k+1));
               Q_mixed_derivs.set(4,hx+1,hy+1,k+2,0.0);
               Q_mixed_derivs.set(5,hx+1,hy+1,k+2,-dt_dx*((double)hx+1.0)/(k+1)*Ge.get(1,hx+2,hy+1,k+1)-dt_dy*((double)hy+1.0)/(k+1)*Ge.get(2,hx+1,hy+2,k+1));
           }

       }

        for(int k=MAX_FLUX_DERIVS-1;k<MAX_DERIVS;k++)
        {

          //Populate time derivatives of (rho u)^2...
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                double tmp3 = 0.0;
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*Q_mixed_derivs.get(2, hx+1-rx,hy+1-ry, k+1-s );
                      tmp2 += Q_mixed_derivs.get( 3, rx+1,ry+1, s+1)*Q_mixed_derivs.get(3, hx+1-rx,hy+1-ry, k+1-s );
                      tmp3 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*Q_mixed_derivs.get(3, hx+1-rx,hy+1-ry, k+1-s );
                }
                Gs.set( 1, hx+1,hy+1, k+1, tmp1 );
                Gs.set( 2, hx+1,hy+1, k+1, tmp2 );
                Gs.set( 3, hx+1,hy+1, k+1, tmp3 );
          }

          //Populate time derivatives of 1/rho...
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                //if(hx+hy>0)
                {

                double tmp = 0.0;
                R.set(1,hx+1,hy+1,k+1,0.0);
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1,1));
                R.set( 1, hx+1,hy+1, k+1, tmp );
                }
          }

          //populate the time derivatives of (rho u)^2/rho and u
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp1  = 0.0;
                double tmp2  = 0.0;
                double tmp3  = 0.0;
                double tmp11 = 0.0;
                double tmp12 = 0.0;
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += Gs.get( 1, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp2 += Gs.get( 2, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp3 += Gs.get( 3, rx+1,ry+1, s+1)*R.get(1, hx+1-rx, hy+1-ry,k+1-s );
                      tmp11 += Q_mixed_derivs.get( 2, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                      tmp12 += Q_mixed_derivs.get( 3, rx+1,ry+1, s+1)*R.get(1, hx+1-rx,hy+1-ry, k+1-s );
                }
                G2r.set( 1, hx+1,hy+1, k+1, tmp1 );
                G2r.set( 2, hx+1,hy+1, k+1, tmp2 );
                G2r.set( 3, hx+1,hy+1, k+1, tmp3 );
                G1r.set( 1, hx+1,hy+1, k+1, tmp11 );
                G1r.set( 2, hx+1,hy+1, k+1, tmp12 );
          }

          //Populate time derivatives of pressure
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp = 0.0;
                tmp = (gamma-1.0)*(Q_mixed_derivs.get( 5, hx+1,hy+1, k+1)-0.5*G2r.get(1,hx+1,hy+1,k+1)-0.5*G2r.get(2,hx+1,hy+1,k+1));
                P.set( 1, hx+1,hy+1, k+1, tmp );
          }
       
          //populate time derivatives of energy equation flux...
          for( int hx=0; hx < MAX_DERIVS; hx++ )
          for( int hy=0; hy < MAX_DERIVS; hy++ )
          {
                double tmp1 = 0.0;
                double tmp2 = 0.0;
                for( int rx=0; rx <= hx; rx++ )
                for( int ry=0; ry <= hy; ry++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp1 += (Q_mixed_derivs.get( 5, rx+1,ry+1, s+1)+P.get( 1, rx+1,ry+1, s+1))*G1r.get(1, hx+1-rx,hy+1-ry, k+1-s );
                      tmp2 += (Q_mixed_derivs.get( 5, rx+1,ry+1, s+1)+P.get( 1, rx+1,ry+1, s+1))*G1r.get(2, hx+1-rx,hy+1-ry, k+1-s );
                }
                Ge.set(1,hx+1,hy+1,k+1,tmp1);
                Ge.set(2,hx+1,hy+1,k+1,tmp2);

           }

        }
       
        // ------------------------------------------------------------------ //
        //
        // Construct the time-averaged flux.
        //
        // Now, construct time-averaged flux by computing a numerical
        // approximation to the integral of the flux function, given a
        // space-time expansion of the conserved variables.  By computing this
        // term, the scheme stays automatically conservative.
        //
        // Note that Q_mixed_derivs already contains a factorial in its 
        // definition.  In order to
        // construct the "time-averaged" version, we need to only divide by a
        // single extra factor of (k+1).  
        //
        // ------------------------------------------------------------------ //

        // 1nd-component of flux function : F in q_t + F_x + G_y = Psi
        double rho_ta_flux  = 0.;
        double u1_ta_flux   = 0.;
        double u2_ta_flux   = 0.;
        double E_ta_flux    = 0.;
        for( int mq=1; mq <= NumQuad; mq++ )
        {
            // Evaluate q at this quadrature point
            double rho  = Q_mixed_derivs.get(1,1,1,1);
            double ru1  = Q_mixed_derivs.get(2,1,1,1);
            double ru2  = Q_mixed_derivs.get(3,1,1,1);
            double E    = Q_mixed_derivs.get(5,1,1,1);
            double P1   = P.get(1,1,1,1);
            double Ge1  = Ge.get(1,1,1,1);
            for( int k=1; k < MAX_DERIVS; k++ )
            { 
                rho += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 1, 1,1, k+1 ); 
                ru1 += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 2, 1,1, k+1 );
                ru2 += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 3, 1,1, k+1 ); 
                E   += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 5, 1,1, k+1 );

                //the pressure function 
                P1 += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*P.get( 1, 1,1, k+1 ); 

                //the function that gives the flux for the energy equation
                Ge1 += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Ge.get( 1, 1,1, k+1 ); 
                //{printf(" here rho=%le, rhou1=%le, rhou2=%le, P1=%le,E=%le \n",rho,ru1,ru2,P1,E);}

            }

            //printf("here F1=%le, F2=%le, F3=%le, F5=%le \n",rho_ta_flux,u1_ta_flux,u2_ta_flux,E_ta_flux);
            //  printf("THIS %le \n",w1d.get( mq ));
            rho_ta_flux += (0.5*w1d.get( mq )) * ( ru1 );
            u1_ta_flux  += (0.5*w1d.get( mq )) * ( ru1*ru1/rho+P1);
            u2_ta_flux  += (0.5*w1d.get( mq )) * ( ru1*ru2/rho);
            E_ta_flux   += (0.5*w1d.get( mq )) * (E+P1)*ru1/rho;
            //printf("here F1=%le, F2=%le, F3=%le, F5=%le \n",rho_ta_flux,u1_ta_flux,u2_ta_flux,E_ta_flux);

        }

        F.set( i,j, 1, rho_ta_flux ); // flux for density
        F.set( i,j, 2, u1_ta_flux );  // flux for 1-component of momentum
        F.set( i,j, 3, u2_ta_flux );  // flux for 2-component of momentum
        F.set( i,j, 4, 0.0 );         // not used in 1/2D (flux for 3-component of momentum)
        F.set( i,j, 5, E_ta_flux );   // flux for energy


        // 2nd-component of flux function : G in q_t + F_x + G_y = Psi
        rho_ta_flux  = 0.;
        u1_ta_flux   = 0.;
        u2_ta_flux   = 0.;
        E_ta_flux    = 0.;
        for( int mq=1; mq <= NumQuad; mq++ )
        {

            // Evaluate q at this quadrature point
            double rho  = Q_mixed_derivs.get(1,1,1,1);
            double ru1  = Q_mixed_derivs.get(2,1,1,1);
            double ru2  = Q_mixed_derivs.get(3,1,1,1);
            double E    = Q_mixed_derivs.get(5,1,1,1);
            double P1   = P.get(1,1,1,1);
            double Ge1  = Ge.get(2,1,1,1);
            for( int k=1; k < MAX_DERIVS; k++ )
            { 
                rho += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 1, 1,1, k+1 ); 
                ru1 += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 2, 1,1, k+1 );
                ru2 += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 3, 1,1, k+1 ); 
                E   += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 5, 1,1, k+1 );

                //the pressure function 
                P1 += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*P.get( 1, 1,1, k+1 ); 

                //the function that gives the flux for the energy equation
                Ge1 += ( pow(0.5*( 1.0 + x1d.get(mq) ), k) )*Ge.get( 2, 1,1, k+1 ); 
            }

            // {printf(" here %le %le \n",rho,P1);}
            rho_ta_flux += (0.5*w1d.get( mq )) * ( ru2 );
            u1_ta_flux  += (0.5*w1d.get( mq )) * ( ru2*ru1/rho);
            u2_ta_flux  += (0.5*w1d.get( mq )) * ( ru2*ru2/rho+P1);
            E_ta_flux   += (0.5*w1d.get( mq )) * (E+P1)*ru2/rho;
            // printf("here G1=%le, G2=%le, G3=%le, G5=%le \n",rho_ta_flux,u1_ta_flux,u2_ta_flux,E_ta_flux);

        }

        G.set( i,j, 1, rho_ta_flux ); // flux for density
        G.set( i,j, 2, u1_ta_flux );  // flux for 1-component of momentum
        G.set( i,j, 3, u2_ta_flux );  // flux for 2-component of momentum
        G.set( i,j, 4, 0.0 );         // not used in 1/2D (flux for 3-component of momentum)
        G.set( i,j, 5, E_ta_flux );   // flux for energy

    }

}
