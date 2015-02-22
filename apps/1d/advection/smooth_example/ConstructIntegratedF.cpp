#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"
#include "StateVars.h"
#include "assert.h"

// Time-integrated flux-function.  This one uses new "Diff1" and "Diff2" to try
// and define a PIF-WENO method that is TVD below a certain threshold on the CFL
// number.
//
// See: "The Picard integral formulation of weighted essentially non-oscillatory
// schemes" available on arxiv: http://arxiv.org/abs/1403.1282

// Central Finite difference approximations.  See lib/WenoReconstruct.cpp
void Diff1      ( double dx, const dTensor2& f, dTensor1& fx );
void Diff2      ( double dx, const dTensor2& f, dTensor1& fxx );
//void Diff1NC    ( double dx, const dTensor2& f, dTensor1& fx );
//void Diff2NC    ( double dx, const dTensor2& f, dTensor1& fxx );

void FluxFunc  (const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
void DFluxFunc(const dTensor1& xpts, const dTensor2& Q, const dTensor2& Aux,
    dTensor3& Dflux);
void D2FluxFunc(const dTensor1& xpts, const dTensor2& Q, const dTensor2& Aux,
    dTensor4& D2flux);

// Used for construcing the flux function
void SampleFunction( 
    int istart, int iend,
    const dTensorBC2& qin, 
    const dTensorBC2& auxin,  
    dTensorBC2& Fout,
    void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

// This function computes the (linear) finite difference approximation to the
// integrated flux function on the conserved variables.  It requires knowledge
// of FluxFunc (1st-order), DFluxFunc (2nd-order) and D2FluxFunc (3rd-order)
// in order to compute the expansion:
//
//     F := f - dt/2 A f_x 
//            + dt^2/3 ( A_q ( f_x, f_x ) + A ( A_q ( q_x, f_x ) + A f_xx )
//            + \cdots.
//
// Where the flux Jacobian and Hessian are defined as:
//
//      A := \partial f   / \partial q,  and 
//    A_q := \partial^2 f / \partial^2 q.
//
// Higher order methods would require further terms to be defined here,
// including "super"-Hessians.
//
// Lax-Friedrichs flux splitting + WENO reconstruction can then be applied 
// to the integrated flux function, F to define an update of the form:
//
//     q^{n+1} = q^n - \dt F_x.
//
// See also: DFluxFunc and D2FluxFunc.
void ConstructIntegratedF( double dt, 
    const StateVars& Q,
    dTensorBC1& smax, dTensorBC2& F)
{

    const dTensorBC2& q   = Q.const_ref_q();
    const dTensorBC2& aux = Q.const_ref_aux();

// TODO - this shouldn't really be placed here ... -DS
const double speed = aux.get(1,1);
const double s2    = speed*speed;
const double s3    = speed*s2;

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    // Compute the first derivative, q_x (and/or f_x)
    dTensorBC2 q_x (mx, meqn, mbc );
    dTensorBC2 q_xx(mx, meqn, mbc );

// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );

    // Compute finite difference approximation to q_x
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

        // First derivative
        dTensor1 qx_val  ( meqn );
        // Diff1NC( dx, qvals, qx_val  );
        Diff1( dx, qvals, qx_val  );

        // Compute the product: A f_x
        for( int m1=1; m1 <= meqn; m1++ )
        {
            q_x.set( i, m1, qx_val.get( m1 ) );
        }

    }


    // Compute the second-derivative (using non-conservative WENO
    // differentiation)
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    {

        // Physical location for this current value:
        dTensor1 xpts( 1 );
        xpts.set( 1, xlow + double(i)*dx - 0.5*dx );

        // Second derivative
        dTensor2 qx_vals( meqn, mpts_sten );
        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                qx_vals.set( m, r, q_x.get( i+s, m ) );
                s++;
            }
        }

        // Second derivative
        dTensor1 qxx_val  ( meqn );
        // Diff1NC( dx, qx_vals, qxx_val  );
        Diff1( dx, qx_vals, qxx_val  );

        for( int m=1; m<=meqn; m++ )
        {
            F.set( i, m, speed*q.get(i,m) - 0.5*dt*(s2*q_x.get(i,m) - dt/3.*s3*qxx_val.get(m) ) );
        }

    }

}

void LocalIntegrate( 
    int nterms,
    double dx, double xc, int meqn, int maux, int mpts_sten, int half_mpts_sten,
    const int i, const dTensorBC2& q, const dTensorBC2& aux, const dTensorBC2& f, dTensor1& f_t, dTensor1& f_tt )
{

        // Physical location for this current value:
        dTensor1 xpts( 1 );
        xpts.set( 1, xc );

        // Save the flux function:
        dTensor2 Fvals( meqn, mpts_sten );
        dTensor2 qvals( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                Fvals.set( m, r, f.get( i+s, m ) );
                qvals.set( m, r, q.get( i+s, m ) );
                s++;
            }
        }

        // Storage for local derivatives:
        dTensor1 fx_val  ( meqn );
        dTensor1 qx_val  ( meqn );
        dTensor1 fxx_val ( meqn );

        // Compute a FD approximation to the derivatives:
        Diff1( dx, Fvals, fx_val  );
        Diff1( dx, qvals, qx_val  );
        Diff2( dx, Fvals, fxx_val );

//      Diff1NC( dx, Fvals, fx_val  );
//      Diff1NC( dx, qvals, qx_val  );
//      Diff2NC( dx, Fvals, fxx_val );

        // Construct the first product: A f_x:
        dTensor3 A( 1, meqn, meqn );
        dTensor2 q_transpose( 1, meqn );
        dTensor2 a_transpose( 1, maux );

        for( int m=1; m <= meqn; m++ )
        {
            q_transpose.set(1, m, q.get(i,m) );
        }
        for( int m=1; m <= maux; m++ )
        {
            a_transpose.set(1, m, aux.get(i,m) );
        }

        // Compute the Jacobian:
        DFluxFunc(xpts, q_transpose, a_transpose, A);

        // Compute the product: A f_x
        for( int m1=1; m1 <= meqn; m1++ )
        {
            double tmp = 0.;
            for( int m2=1; m2 <= meqn; m2++ )
            {
                tmp += -A.get(1, m1, m2) * fx_val.get(m2);
            }
            f_t.set( m1, tmp );
        }

        // ---  Third-order terms --- //
        if( nterms > 2 )
        {

            // Hessian
            dTensor4 H( 1, meqn, meqn, meqn );
            D2FluxFunc(xpts, q_transpose, a_transpose, H);

            // Compute terms that get multiplied by \pd2{ f }{ q }.
            dTensor1 tmp_vec( meqn );
            for( int m =1; m <= meqn; m++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;
                for( int m1=1; m1 <= meqn; m1++ )
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += H.get(1,m,m1,m2) * fx_val.get(m1)*fx_val.get(m2);
                    tmp2 += H.get(1,m,m1,m2) * qx_val.get(m1)*fx_val.get(m2);
                }
                f_tt.set( m, tmp1 );
                tmp_vec.set( m, tmp2 );
            }

            // Add in the third term that gets multiplied by A:
            for( int m1=1; m1 <= meqn; m1++ )
            {
                double tmp = 0.;
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp += A.get(1,m1,m2)*fxx_val.get(m2);
                }
                tmp_vec.set( m1, tmp_vec.get(m1) + tmp );
            }

            // Multiply final term by A:
            for( int m1=1; m1 <= meqn; m1++ )
            {
                double tmp = 0.;
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp += A.get(1,m1,m2)*tmp_vec.get(m2);
                }
                f_tt.set( m1, f_tt.get(m1) + tmp );
            }

        }

}

void ConstructIntegratedF( double dt, 
    double alpha1, double beta1, const StateVars& Q1,
    double alpha2, double beta2, const StateVars& Q2,
    dTensorBC1& smax, dTensorBC2& F)
{

    const dTensorBC2& q1 = Q1.const_ref_q();
    const dTensorBC2& q2 = Q2.const_ref_q();

    const dTensorBC2& aux1 = Q1.const_ref_aux();
    const dTensorBC2& aux2 = Q2.const_ref_aux();

    const int mx     = q1.getsize(1);
    const int meqn   = q1.getsize(2);
    const int maux   = aux1.getsize(2);
    const int mbc    = q1.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC2 f1( mx, meqn, mbc );
    dTensorBC2 f2( mx, meqn, mbc );
    SampleFunction( 1-mbc, mx+mbc, q1, aux1, f1, &FluxFunc );
    SampleFunction( 1-mbc, mx+mbc, q2, aux2, f2, &FluxFunc );

// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    {

        dTensor1 f1_t( meqn );    f1_t.setall(0.);
        dTensor1 f1_tt( meqn );   f1_tt.setall(0.);
        dTensor1 f2_t( meqn );    f2_t.setall( 0.);
        dTensor1 f2_tt( meqn );   f2_tt.setall(0.);

        double xc = xlow + double(i)*dx - 0.5*dx;
        LocalIntegrate( 2, dx, xc, meqn, maux, mpts_sten, half_mpts_sten,
            i, q1, aux1, f1, f1_t, f1_tt );

        LocalIntegrate( 2, dx, xc, meqn, maux, mpts_sten, half_mpts_sten,
            i, q2, aux2, f2, f2_t, f2_tt );

        // Second/Third-order accuracy:
        for( int m=1; m<=meqn; m++ )
        {
            F.set( i, m, alpha1*f1.get(i,m) + beta1*dt*(f1_t.get(m)) + 
                         alpha2*f2.get(i,m) + beta2*dt*(f2_t.get(m)) );
        }


    }

}

void ConstructIntegratedF( double dt, 
    double alpha1, double beta1, const StateVars& Q1,
    double alpha2, double beta2, const StateVars& Q2,
    double alpha3, double beta3, const StateVars& Q3,
    dTensorBC1& smax, dTensorBC2& F)
{

    const dTensorBC2& q1 = Q1.const_ref_q();
    const dTensorBC2& q2 = Q2.const_ref_q();
    const dTensorBC2& q3 = Q3.const_ref_q();

    const dTensorBC2& aux1 = Q1.const_ref_aux();
    const dTensorBC2& aux2 = Q2.const_ref_aux();
    const dTensorBC2& aux3 = Q3.const_ref_aux();



    const int mx     = q1.getsize(1);
    const int meqn   = q1.getsize(2);
    const int maux   = aux1.getsize(2);
    const int mbc    = q1.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC2 f1( mx, meqn, mbc );
    dTensorBC2 f2( mx, meqn, mbc );
    dTensorBC2 f3( mx, meqn, mbc );
    SampleFunction( 1-mbc, mx+mbc, q1, aux1, f1, &FluxFunc );
    SampleFunction( 1-mbc, mx+mbc, q2, aux2, f2, &FluxFunc );
    SampleFunction( 1-mbc, mx+mbc, q3, aux3, f3, &FluxFunc );

// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    {

        dTensor1 f1_t( meqn );    f1_t.setall(0.);
        dTensor1 f1_tt( meqn );   f1_tt.setall(0.);
        dTensor1 f2_t( meqn );    f2_t.setall( 0.);
        dTensor1 f2_tt( meqn );   f2_tt.setall(0.);
        dTensor1 f3_t( meqn );    f3_t.setall( 0.);
        dTensor1 f3_tt( meqn );   f3_tt.setall(0.);

        double xc = xlow + double(i)*dx - 0.5*dx;
        LocalIntegrate( 2, dx, xc, meqn, maux, mpts_sten, half_mpts_sten,
            i, q1, aux1, f1, f1_t, f1_tt );

        LocalIntegrate( 2, dx, xc, meqn, maux, mpts_sten, half_mpts_sten,
            i, q2, aux2, f2, f2_t, f2_tt );

        LocalIntegrate( 2, dx, xc, meqn, maux, mpts_sten, half_mpts_sten,
            i, q3, aux3, f3, f3_t, f3_tt );

        // Second/Third-order accuracy:
        for( int m=1; m<=meqn; m++ )
        {
            F.set( i, m, alpha1*f1.get(i,m) + beta1*dt*(f1_t.get(m)) + 
                         alpha2*f2.get(i,m) + beta2*dt*(f2_t.get(m)) +
                         alpha3*f3.get(i,m) + beta3*dt*(f3_t.get(m)) );
        }

    }

}

// Two-stage, three derivative method
void ConstructIntegratedF( double dt, 
    double alpha1, double beta1, double charlie1, const StateVars& Q1,
    double alpha2, double beta2, double charlie2, const StateVars& Q2,
    dTensorBC1& smax, dTensorBC2& F)
{

    const dTensorBC2& q1 = Q1.const_ref_q();
    const dTensorBC2& q2 = Q2.const_ref_q();

    const dTensorBC2& aux1 = Q1.const_ref_aux();
    const dTensorBC2& aux2 = Q2.const_ref_aux();


    const int mx     = q1.getsize(1);
    const int meqn   = q1.getsize(2);
    const int maux   = aux1.getsize(2);
    const int mbc    = q1.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC2 f1( mx, meqn, mbc );
    dTensorBC2 f2( mx, meqn, mbc );
    SampleFunction( 1-mbc, mx+mbc, q1, aux1, f1, &FluxFunc );
    SampleFunction( 1-mbc, mx+mbc, q2, aux2, f2, &FluxFunc );

// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    {

        dTensor1 f1_t( meqn );    f1_t.setall(0.);
        dTensor1 f1_tt( meqn );   f1_tt.setall(0.);
        dTensor1 f2_t( meqn );    f2_t.setall( 0.);
        dTensor1 f2_tt( meqn );   f2_tt.setall(0.);

        double xc = xlow + double(i)*dx - 0.5*dx;
        LocalIntegrate( 3, dx, xc, meqn, maux, mpts_sten, half_mpts_sten,
            i, q1, aux1, f1, f1_t, f1_tt );

        LocalIntegrate( 3, dx, xc, meqn, maux, mpts_sten, half_mpts_sten,
            i, q2, aux2, f2, f2_t, f2_tt );

        // Second/Third-order accuracy:
        for( int m=1; m<=meqn; m++ )
        {
            F.set( i, m, 
                alpha1*f1.get(i,m) + beta1*dt*(f1_t.get(m)) + charlie1*dt*dt*f1_tt.get(m) +
                alpha2*f2.get(i,m) + beta2*dt*(f2_t.get(m)) + charlie2*dt*dt*f2_tt.get(m) );
        }


    }

}
