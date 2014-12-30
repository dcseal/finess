#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "assert.h"
#include "IniParams.h"
#include "StateVars.h"

// Central difference formulae
// see $FINESS/lib/WenoReconstruct.cpp
void Diff1( double dx, const dTensor2& f, dTensor1& fx );
void Diff2( double dx, const dTensor2& f, dTensor1& fxx );
void SampleFunctionTypeB( 
        int istart, int iend,
        int jstart, int jend,
        int kstart, int kend,
        const dTensorBC4& qin, 
        const dTensorBC4& auxin,  
        dTensorBC5& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&));

void FluxFunc(const dTensor2& xpts, const dTensor2& Q,const dTensor2& Aux, dTensor3& flux);

void DFluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux,
    dTensor4& Dflux);
void D2FluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux,
    dTensor5& D2flux);

void ConstructIntegratedR( double dt, const StateVars& Q,
        dTensorBC4& smax, dTensorBC4& F, dTensorBC4& G, dTensorBC4& H)
{

    const dTensorBC4& q   = Q.const_ref_q  ();
    const dTensorBC4& aux = Q.const_ref_aux();

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int mz     = global_ini_params.get_mz();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double dy    = global_ini_params.get_dy();
    const double dz    = global_ini_params.get_dz();
    const double xlow  = global_ini_params.get_xlow();
    const double ylow  = global_ini_params.get_ylow();
    const double zlow  = global_ini_params.get_zlow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC5 R( mx, my, mz, meqn, 3, mbc );  // place-holder for the flux function
    SampleFunctionTypeB( 1-mbc, mx+mbc, 1-mbc, my+mbc, 1-mbc, mz+mbc, q, aux, R, &FluxFunc );

    // TODO  - allow for different sized stencils for different orders (-DS)
    const int mbc_small      = 3;
    const int      mpts_sten = 5;
    //const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );
    const int half_mpts_sten = 3;

    const int ndim = 3;

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
        for( int j = 1-mbc_small; j <= my+mbc_small; j++ )
            for( int k = 1-mbc_small; k <= mz+mbc_small; k++ )
            {

                // Physical location for this current value:
                dTensor2 xpts( 1, ndim );
                xpts.set( 1, 1, xlow + double(i)*dx - 0.5*dx );
                xpts.set( 1, 2, ylow + double(j)*dy - 0.5*dy );
                xpts.set( 1, 3, zlow + double(k)*dz - 0.5*dz );

                // Save the flux function:
                dTensor2 Fvals  ( meqn, mpts_sten );
                dTensor2 Gvals  ( meqn, mpts_sten );
                dTensor2 Hvals  ( meqn, mpts_sten );
                dTensor2 qvalsx ( meqn, mpts_sten );
                dTensor2 qvalsy ( meqn, mpts_sten );
                dTensor2 qvalsz ( meqn, mpts_sten );

                for( int m=1; m <= meqn; m++ )
                {
                    int s = -half_mpts_sten+1;
                    for( int r = 1; r <= mpts_sten; r++ )
                    {
                        Fvals.set( m, r, R.get( i+s, j, k, m, 1 ) ); // for computing f_x
                        Gvals.set( m, r, R.get( i, j+s, k, m, 2 ) ); // for computing g_y
                        Hvals.set( m, r, R.get( i, j, k+s, m, 3 ) ); // for computing h_z
                        qvalsx.set( m, r, q.get( i+s, j, k, m ) );
                        qvalsy.set( m, r, q.get( i, j+s, k, m ) );
                        qvalsz.set( m, r, q.get( i, j, k+s, m ) );
                        s++;
                    }
                }

                // Storage for local derivatives:
                dTensor1 fx_val  ( meqn ), qx_val( meqn );
                dTensor1 gy_val  ( meqn ), qy_val( meqn );
                dTensor1 hz_val  ( meqn ), qz_val( meqn );

                // Compute a FD approximation to the derivatives:
                Diff1( dx, Fvals,  fx_val  );
                Diff1( dx, qvalsx, qx_val  );

                Diff1( dy, Gvals,  gy_val  );
                Diff1( dy, qvalsy, qy_val  );

                Diff1( dz, Hvals,  hz_val  );
                Diff1( dz, qvalsz, qz_val  );

                // A(1, :, :, 1)  -- f'(q)
                // A(1, :, :, 2)  -- g'(q)
                // A(1, :, :, 3)  -- h'(q)
                dTensor4 A( 1, meqn, meqn, 3 ); 
                dTensor2 q_transpose( 1, meqn );
                dTensor2 a_transpose( 1, maux );

                for( int m=1; m <= meqn; m++ )
                {
                    q_transpose.set(1, m, q.get(i,j,k,m) );
                }
                for( int m=1; m <= maux; m++ )
                {
                    a_transpose.set(1, m, aux.get(i,j,k,m) );
                }

                // Compute the Jacobian:
                DFluxFunc(xpts, q_transpose, a_transpose, A);

                // Compute the product: f'(q)*(f_x+g_y) + g'(q)*(f_x+g_y)
                dTensor1 f_t( meqn ), g_t( meqn ), h_t(meqn);
                for( int m1=1; m1 <= meqn; m1++ )
                {
                    double tmp1 = 0.;
                    double tmp2 = 0.;
                    double tmp3 = 0.;
                    for( int m2=1; m2 <= meqn; m2++ )
                    {
                        tmp1 += -(A.get(1, m1, m2, 1)) * ( fx_val.get(m2) + gy_val.get(m2) + hz_val.get(m2));
                        tmp2 += -(A.get(1, m1, m2, 2)) * ( fx_val.get(m2) + gy_val.get(m2) + hz_val.get(m2));
                        tmp3 += -(A.get(1, m1, m2, 3)) * ( fx_val.get(m2) + gy_val.get(m2) + hz_val.get(m2));
                    }
                    f_t.set( m1, tmp1 );
                    g_t.set( m1, tmp2 );
                    h_t.set( m1, tmp3 );

                }

                // ---  Third-order terms --- //
                dTensor1 f_tt( meqn );   f_tt.setall(0.);
                dTensor1 g_tt( meqn );   g_tt.setall(0.);
                dTensor1 h_tt( meqn );   h_tt.setall(0.);
                if( global_ini_params.get_time_order() > 2 )
                {

                    // Hessian
                    dTensor5 H( 1, meqn, meqn, meqn, 3 );
                    D2FluxFunc(xpts, q_transpose, a_transpose, H);

                    // ----------------------------------- //
                    // Part I: Compute (f_x + g_y + h_z)_{,t}
                    // ----------------------------------- //

                    dTensor1 fx_plus_gy_plus_hz_t( meqn ); fx_plus_gy_plus_hz_t.setall(0.);

                    // Start with a term that get's used frequently:
                    dTensor1 fx_plus_gy_plus_hz( meqn );
                    for( int m =1; m <= meqn; m++ )
                    {
                        double tmp = fx_val.get(m) + gy_val.get(m) + hz_val.get(m);
                        fx_plus_gy_plus_hz.set(m, tmp );
                    }

                    // Second-derivatives:
                    dTensor1 fxx_val(meqn);
                    dTensor1 fxy_val(meqn);
                    dTensor1 fxz_val(meqn);
                    dTensor1 gxy_val(meqn);
                    dTensor1 gyy_val(meqn);
                    dTensor1 gyz_val(meqn);
                    dTensor1 hxz_val(meqn);
                    dTensor1 hyz_val(meqn);
                    dTensor1 hzz_val(meqn);

                    Diff2( dx, Fvals,  fxx_val  );
                    Diff2( dy, Gvals,  gyy_val  );
                    Diff2( dz, Hvals,  hzz_val  );

                    // Cross - derivaties
                    fxy_val.setall(0.);
                    fxz_val.setall(0.);
                    gxy_val.setall(0.);
                    gyz_val.setall(0.);
                    hxz_val.setall(0.);
                    hyz_val.setall(0.);

                    // 2nd-order stencil (for mixed derivatives)
                    for( int m=1; m <= meqn; m++ )
                    {
                        fxy_val.set(m, 0.25*(R.get(i+1,j+1,k,m,1)-R.get(i-1,j+1,k,m,1)-R.get(i+1,j-1,k,m,1)+R.get(i-1,j-1,k,m,1))/dx/dy);
                        gxy_val.set(m, 0.25*(R.get(i+1,j+1,k,m,2)-R.get(i-1,j+1,k,m,2)-R.get(i+1,j-1,k,m,2)+R.get(i-1,j-1,k,m,2))/dx/dy);
                        fxz_val.set(m, 0.25*(R.get(i+1,j,k+1,m,1)-R.get(i-1,j,k+1,m,1)-R.get(i+1,j,k-1,m,1)+R.get(i-1,j,k-1,m,1))/dx/dz);
                        hxz_val.set(m, 0.25*(R.get(i+1,j,k+1,m,3)-R.get(i-1,j,k+1,m,3)-R.get(i+1,j,k-1,m,3)+R.get(i-1,j,k-1,m,3))/dx/dz);
                        gyz_val.set(m, 0.25*(R.get(i,j+1,k+1,m,2)-R.get(i,j-1,k+1,m,2)-R.get(i,j+1,k-1,m,2)+R.get(i,j-1,k-1,m,2))/dy/dz);
                        hyz_val.set(m, 0.25*(R.get(i,j+1,k+1,m,3)-R.get(i,j-1,k+1,m,3)-R.get(i,j+1,k-1,m,3)+R.get(i,j-1,k-1,m,3))/dy/dz);
                    }

                    // Compute (f_x + g_y + h_z)_t
                    for( int m =1; m <= meqn; m++ )
                    {
                        double tmp = 0.;

                        // Terms that get multiplied by Hessians:
                        for( int m1=1; m1 <= meqn; m1++ )
                            for( int m2=1; m2 <= meqn; m2++ )
                            {

                                tmp += H.get(1,m,m1,m2,1)*qx_val.get(m1)*fx_plus_gy_plus_hz.get(m2);
                                tmp += H.get(1,m,m1,m2,2)*qy_val.get(m1)*fx_plus_gy_plus_hz.get(m2);
                                tmp += H.get(1,m,m1,m2,3)*qz_val.get(m1)*fx_plus_gy_plus_hz.get(m2);
                            }

                        // Terms that get multiplied by Jacobians:
                        for( int m1=1; m1 <= meqn; m1++ )
                        {

                            tmp += A.get(1,m,m1,1)*( fxx_val.get(m1)+gxy_val.get(m1)+hxz_val.get(m1) );
                            tmp += A.get(1,m,m1,2)*( fxy_val.get(m1)+gyy_val.get(m1)+hyz_val.get(m1) );
                            tmp += A.get(1,m,m1,3)*( fxz_val.get(m1)+gyz_val.get(m1)+hzz_val.get(m1) );                        }

                        fx_plus_gy_plus_hz_t.set( m, tmp );
                    }


                    // ----------------------------------- //
                    // Part II: Compute 
                    //      f'(q) * fx_plus_gy_plus_hz_t and 
                    //      g'(q) * fx_plus_gy_plus_hz_t and
                    //      h'(q) * fx_plus_gy_plus_hz_t
                    // ----------------------------------- //

                    // Add in the third term that gets multiplied by A:
                    for( int m1=1; m1 <= meqn; m1++ )
                    {
                        double tmp1 = 0.;
                        double tmp2 = 0.;
                        double tmp3 = 0.;
                        for( int m2=1; m2 <= meqn; m2++ )
                        {
                            tmp1 += A.get(1,m1,m2,1)*fx_plus_gy_plus_hz_t.get(m2);
                            tmp2 += A.get(1,m1,m2,2)*fx_plus_gy_plus_hz_t.get(m2);
                            tmp3 += A.get(1,m1,m2,3)*fx_plus_gy_plus_hz_t.get(m2);
                        }
                        f_tt.set( m1, tmp1 );
                        g_tt.set( m1, tmp2 );
                        h_tt.set( m1, tmp3 );
                    }

                    // ----------------------------------------------- //
                    // Part III: Add in contributions from
                    //      f''(q) * (fx_plus_gy_plus_hz, fx_plus_gy_plus_hz ) and 
                    //      g''(q) * (fx_plus_gy_plus_hz, fx_plus_gy_plus_hz ) and
                    //      h''(q) * (fx_plus_gy_plus_hz, fx_plus_gy_plus_hz ).
                    // ----------------------------------------------- //
                    for( int m =1; m <= meqn; m++ )
                    {
                        double tmp1 = 0.;
                        double tmp2 = 0.;
                        double tmp3 = 0.;
                        // Terms that get multiplied by the Hessian:
                        for( int m1=1; m1 <= meqn; m1++ )
                            for( int m2=1; m2 <= meqn; m2++ )
                            {
                                tmp1 += H.get(1,m,m1,m2,1)*fx_plus_gy_plus_hz.get(m1)*fx_plus_gy_plus_hz.get(m2);
                                tmp2 += H.get(1,m,m1,m2,2)*fx_plus_gy_plus_hz.get(m1)*fx_plus_gy_plus_hz.get(m2);
                                tmp3 += H.get(1,m,m1,m2,3)*fx_plus_gy_plus_hz.get(m1)*fx_plus_gy_plus_hz.get(m2);
                            }

                        f_tt.set( m, f_tt.get(m) + tmp1 );
                        g_tt.set( m, g_tt.get(m) + tmp2 );
                        h_tt.set( m, h_tt.get(m) + tmp3 );
                    }

                }

                // FINAL STEP: save the time-integrated values:
                for( int m=1; m<=meqn; m++ )

                {
                    F.set(i,j,k,m, R.get(i,j,k,m,1) + 0.5*dt*(f_t.get(m) + dt/3.0*f_tt.get(m)) );
                    G.set(i,j,k,m, R.get(i,j,k,m,2) + 0.5*dt*(g_t.get(m) + dt/3.0*g_tt.get(m)) );
                    H.set(i,j,k,m, R.get(i,j,k,m,3) + 0.5*dt*(h_t.get(m) + dt/3.0*h_tt.get(m)) );
                }


            }

}
