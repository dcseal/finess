#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "assert.h"

// EXPERIMENTAL CODE

// See $FINESS/lib/WenoReconstruct.cpp for these central finite difference methods.
void Diff1( double dx, const dTensor2& f, dTensor1& fx );
void Diff2( double dx, const dTensor2& f, dTensor1& fxx );

// This function computes the (linear) finite difference approximation to the
// integrated flux function on the conserved variables.  It requires knowledge
// of FluxFunc (1st-order), DFluxFunc (2nd-order) and D2FluxFunc (3rd-order)
// in order to compute the expansion:
//
//     F := f - dt/2 ( A (f_x + g_y ) )
//            + \cdots.
//
//     G := g - dt/2 ( B (f_x + g_y ) )
//            + \cdots.
//
// TODO - include the extra terms in here.  For now, I will work on the
//        2nd-order method only.  (-DS).
//
// Where the flux Jacobian and Hessian are defined as:
//
//      A := \partial   f / \partial q,   and 
//      B := \partial   g / \partial q,   and 
//    A_q := \partial^2 f / \partial^2 q, and
//    B_q := \partial^2 g / \partial^2 q.
//
// Higher order methods would require further terms to be defined here,
// including "super"-Hessians.
//
// Lax-Friedrichs flux splitting + WENO reconstruction can then be applied 
// to the integrated flux function, F and G to define an update of the form:
//
//     q^{n+1} = q^n - \dt ( F_x + G_y )
//
// See also: DFluxFunc and D2FluxFunc.
void ConstructIntegratedR( double dt, 
    dTensorBC3& aux, dTensorBC3& q,
    dTensorBC3& smax, 
    dTensorBC3& F, dTensorBC3& G)
{


    void FluxFunc(const dTensor2&,const dTensor2&,const dTensor2&,dTensor3&);
    void DFluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux,
        dTensor4& Dflux);
    void D2FluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux,
    	dTensor5& D2flux);

    // Used for construcing the flux function
    void SampleFunction( 
        int istart, int iend,
        int jstart, int jend,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  
              dTensorBC4& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&));

    // Problem dimensions TODO - the boundary data either: 
    //     a) needs one more point, or 
    //     b) needs to double the number of ghost cells, or
    //     c) One-sided differences at the boundary.
    //
    // I would prefer to run with option (c).  (-DS).

    const int mx     = dogParamsCart2.get_mx();
    const int my     = dogParamsCart2.get_my();
    const int meqn   = dogParams.get_meqn();
    const int maux   = dogParams.get_maux();
    const int mbc    = dogParamsCart2.get_mbc();

    // Needed to define derivatives
    const double dx    = dogParamsCart2.get_dx();
    const double dy    = dogParamsCart2.get_dy();
    const double xlow  = dogParamsCart2.get_xlow();
    const double ylow  = dogParamsCart2.get_ylow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC4 R( mx, my, meqn, 2, mbc );  // place-holder for the flux function
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, R, &FluxFunc );

// TODO  - allow for different sized stencils
const int      mpts_sten = 2*mbc-1;  assert_eq( mpts_sten,      5 );
const int half_mpts_sten =   mbc;    assert_eq( half_mpts_sten, 3 );

const int ndim = 2;

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1; i <= mx; i++ )
    for( int j = 1; j <= my; j++ )
    {

        // Physical location for this current value:
        dTensor2 xpts( 1, ndim );
        xpts.set( 1, 1, xlow + double(i)*dx - 0.5*dx );
        xpts.set( 1, 2, ylow + double(j)*dy - 0.5*dy );

        // Save the flux function:
        dTensor2 Fvals  ( meqn, mpts_sten );
        dTensor2 Gvals  ( meqn, mpts_sten );
        dTensor2 qvalsx ( meqn, mpts_sten );
        dTensor2 qvalsy ( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                Fvals.set( m, r, R.get( i+s, j, m, 1 ) ); // for computing f_x
                Gvals.set( m, r, R.get( i, j+s, m, 2 ) ); // for computing g_y
                qvalsx.set( m, r, q.get( i+s, j, m ) );
                qvalsy.set( m, r, q.get( i, j+s, m ) );
                s++;
            }
        }

        // Storage for local derivatives:
        dTensor1 fx_val  ( meqn ), qx_val( meqn );
        dTensor1 gy_val  ( meqn ), qy_val( meqn );

        // Compute a FD approximation to the derivatives:
        Diff1( dx, Fvals,  fx_val  );
        Diff1( dx, qvalsx, qx_val  );

        Diff1( dy, Gvals,  gy_val  );
        Diff1( dy, qvalsy, qy_val  );

        // Construct the first product: A f_x:
        dTensor4 A( 1, meqn, meqn, 2 );             // ( f'(q), g'(q) )
        dTensor2 q_transpose( 1, meqn );
        dTensor2 a_transpose( 1, maux );

        for( int m=1; m <= meqn; m++ )
        {
            q_transpose.set(1, m, q.get(i,j,m) );
        }
        for( int m=1; m <= maux; m++ )
        {
            a_transpose.set(1, m, aux.get(i,j,m) );
        }

        // Compute the Jacobian:
        DFluxFunc(xpts, q_transpose, a_transpose, A);

        // Compute the product: f'(q)*(f_x+g_y) + g'(q)*(f_x+g_y)
        dTensor1 f_t( meqn ), g_t( meqn );
        for( int m1=1; m1 <= meqn; m1++ )
        {
            double tmp1 = 0.;
            double tmp2 = 0.;
            for( int m2=1; m2 <= meqn; m2++ )
            {
                tmp1 += -(A.get(1, m1, m2, 1)) * ( fx_val.get(m2) + gy_val.get(m2));
                tmp2 += -(A.get(1, m1, m2, 2)) * ( fx_val.get(m2) + gy_val.get(m2));
            }
            f_t.set( m1, tmp1 );
            g_t.set( m1, tmp2 );
        }

        // ---  Third-order terms --- //
        dTensor1 f_tt( meqn );   f_tt.setall(0.);
        dTensor1 g_tt( meqn );   g_tt.setall(0.);
        if( dogParams.get_time_order() > 2 )
        {

            // Hessian
            dTensor5 H( 1, meqn, meqn, meqn, 2 );
            D2FluxFunc(xpts, q_transpose, a_transpose, H);

            // ----------------------------------- //
            // Part I: Compute (f_x + g_y)_{,t}
            // ----------------------------------- //

            dTensor1 fx_plus_gy_t( meqn ); fx_plus_gy_t.setall(0.);

            // Start with a term that get's used frequently:
            dTensor1 fx_plus_gy( meqn );
            for( int m =1; m <= meqn; m++ )
            {
                double tmp = fx_val.get(m) + gy_val.get(m);
                fx_plus_gy.set(m, tmp );
            }

            // Second-derivatives:
            dTensor1 fxx_val( meqn ), gyy_val( meqn );
            Diff2( dx, Fvals,  fxx_val  );
            Diff2( dy, Gvals,  gyy_val  );

            // Cross - derivaties
            dTensor1 fxy_val  ( meqn );  fxy_val.setall(0.);
            dTensor1 gxy_val  ( meqn );  gxy_val.setall(0.);
            // TODO - Define the stencil and compute these derivatives!

            // Compute terms that get multiplied by 
            //     \pd2{ f }{ q } and \pd2{ g }{ q }.
            for( int m =1; m <= meqn; m++ )
            {
                double tmp = 0.;

                // Terms that get multiplied by the Hessian:
                for( int m1=1; m1 <= meqn; m1++ )
                for( int m2=1; m2 <= meqn; m2++ )
                {

                    tmp += H.get(1,m,m1,m2,1)*qx_val.get(m1)*fx_plus_gy.get(m2);
                    tmp += H.get(1,m,m1,m2,2)*qy_val.get(m1)*fx_plus_gy.get(m2);
                }

                // Terms that get multiplied by f'(q) and g'(q):
                for( int m1=1; m1 <= meqn; m1++ )
                {

                    tmp += A.get(1,m,m1,1)*( fxx_val.get(m1)+gxy_val.get(m1) );
                    tmp += A.get(1,m,m1,2)*( fxy_val.get(m1)+gyy_val.get(m1) );
                }

                fx_plus_gy_t.set( m, tmp );
            }


            // ----------------------------------- //
            // Part II: Compute 
            //      f'(q) * fx_plus_gy_t and 
            //      g'(q) * fx_plus_gy_t
            // ----------------------------------- //

            // Add in the third term that gets multiplied by A:
            for( int m1=1; m1 <= meqn; m1++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += A.get(1,m1,m2,1)*fx_plus_gy_t.get(m2);
                    tmp2 += A.get(1,m1,m2,2)*fx_plus_gy_t.get(m2);
                }
                f_tt.set( m1, tmp1 );
                g_tt.set( m1, tmp2 );
            }

            // ----------------------------------------------- //
            // Part III: Add in contributions from
            //      f''(q) * (fx_plus_gy, fx_plus_gy ) and 
            //      g''(q) * (fx_plus_gy, fx_plus_gy ).
            // ----------------------------------------------- //
            for( int m =1; m <= meqn; m++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;

                // Terms that get multiplied by the Hessian:
                for( int m1=1; m1 <= meqn; m1++ )
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += H.get(1,m,m1,m2,1)*fx_plus_gy.get(m1)*fx_plus_gy.get(m2);
                    tmp2 += H.get(1,m,m1,m2,2)*fx_plus_gy.get(m1)*fx_plus_gy.get(m2);
                }

                f_tt.set( m, f_tt.get(m) + tmp1 );
                g_tt.set( m, g_tt.get(m) + tmp2 );
            }

        }

        // FINAL STEP: save the time-integrated values:
        for( int m=1; m<=meqn; m++ )
        {
            F.set(i,j,m, R.get(i,j,m,1) + 0.5*dt*(f_t.get(m) + dt/3.0*f_tt.get(m)) );
            G.set(i,j,m, R.get(i,j,m,2) + 0.5*dt*(g_t.get(m) + dt/3.0*g_tt.get(m)) );
        }


    }

    // TODO - something needs to be done about the boundary data!!!
    // For now, we'll assume F satisfies the same boundary conditions that Q
    // does (but this is not true!!)
//void SetBndValues(dTensorBC3&, dTensorBC3&);
//SetBndValues(aux, F );
//SetBndValues(aux, G );

}
