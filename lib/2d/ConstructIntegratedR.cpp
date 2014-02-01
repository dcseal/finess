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
    dTensorBC2& smax, 
    dTensorBC3& F, dTensorBC3& G)
{


    void FluxFunc(const dTensor2&,const dTensor2&,const dTensor2&,dTensor3&);
    void DFluxFunc(const dTensor1& xpts, const dTensor2& Q, const dTensor2& Aux,
	       dTensor3& Dflux);
    void D2FluxFunc(const dTensor1& xpts, const dTensor2& Q, const dTensor2& Aux,
		dTensor4& D2flux);

    // Used for construcing the flux function
    void SampleFunction( 
        int istart, int iend,
        int jstart, int jend,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  
              dTensorBC4& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&));

    // Problem dimensions (TODO - the boundary data either a) needs one more
    // point, or b) needs to double the number of ghost cells )
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
        xpts.set( 1, 2, ylow + double(i)*dy - 0.5*dy );

        // Save the flux function:
        dTensor2 Fvals( meqn, mpts_sten );
        dTensor2 Gvals( meqn, mpts_sten );
        dTensor2 qvals( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                Fvals.set( m, r, R.get( i+s, j, m, 1 ) );
                Gvals.set( m, r, R.get( i+s, j, m, 2 ) );
                qvals.set( m, r, q.get( i+s, j, m ) );
                s++;
            }
        }

        // Storage for local derivatives:
        // TODO - figure out what needs to be stored for the Integrated flux ...
//      dTensor1 fx_val  ( meqn );
//      dTensor1 fy_val  ( meqn );
//      dTensor1 qx_val  ( meqn );
//      dTensor1 qy_val  ( meqn );
//      dTensor1 fxx_val ( meqn );
//      dTensor1 gyy_val ( meqn );

//      // Compute a FD approximation to the derivatives:
//      Diff1( dx, Fvals, fx_val  );
//      Diff1( dx, qvals, qx_val  );
//      Diff2( dx, Fvals, fxx_val );

//      // Construct the first product: A f_x:
//      dTensor3 A( 1, meqn, meqn );
//      dTensor2 q_transpose( 1, meqn );
//      dTensor2 a_transpose( 1, maux );

//      for( int m=1; m <= meqn; m++ )
//      {
//          q_transpose.set(1, m, q.get(i,m) );
//      }
//      for( int m=1; m <= maux; m++ )
//      {
//          a_transpose.set(1, m, aux.get(i,m) );
//      }

//      // Compute the Jacobian:
//      DFluxFunc(xpts, q_transpose, a_transpose, A);

//      // Compute the product: A f_x
//      dTensor1 f_t( meqn );
//      for( int m1=1; m1 <= meqn; m1++ )
//      {
//          double tmp = 0.;
//          for( int m2=1; m2 <= meqn; m2++ )
//          {
//              tmp += -A.get(1, m1, m2) * fx_val.get(m2);
//          }
//          f_t.set( m1, tmp );
//      }

//      // ---  Third-order terms --- //
//      dTensor1 f_tt( meqn );   f_tt.setall(0.);
//      if( dogParams.get_time_order() > 2 )
//      {

//          // Hessian
//          dTensor4 H( 1, meqn, meqn, meqn );
//          D2FluxFunc(xpts, q_transpose, a_transpose, H);

//          // Compute terms that get multiplied by \pd2{ f }{ q }.
//          dTensor1 tmp_vec( meqn );
//          for( int m =1; m <= meqn; m++ )
//          {
//              double tmp1 = 0.;
//              double tmp2 = 0.;
//              for( int m1=1; m1 <= meqn; m1++ )
//              for( int m2=1; m2 <= meqn; m2++ )
//              {
//                  tmp1 += H.get(1,m,m1,m2) * fx_val.get(m1)*fx_val.get(m2);
//                  tmp2 += H.get(1,m,m1,m2) * qx_val.get(m1)*fx_val.get(m2);
//              }
//              f_tt.set( m, tmp1 );
//              tmp_vec.set( m, tmp2 );
//          }

//          // Add in the third term that gets multiplied by A:
//          for( int m1=1; m1 <= meqn; m1++ )
//          {
//              double tmp = 0.;
//              for( int m2=1; m2 <= meqn; m2++ )
//              {
//                  tmp += A.get(1,m1,m2)*fxx_val.get(m2);
//              }
//              tmp_vec.set( m1, tmp_vec.get(m1) + tmp );
//          }

//          // Multiply final term by A:
//          for( int m1=1; m1 <= meqn; m1++ )
//          {
//              double tmp = 0.;
//              for( int m2=1; m2 <= meqn; m2++ )
//              {
//                  tmp += A.get(1,m1,m2)*tmp_vec.get(m2);
//              }
//              f_tt.set( m1, f_tt.get(m1) + tmp );
//          }


//      }

//      // Second/Third-order accuracy:
//      for( int m=1; m<=meqn; m++ )
//      {
//          F.set( i, m, f.get(i,m) + 0.5*dt*(f_t.get(m) + dt/3.0*f_tt.get(m)) );
//      }


    }

    // TODO - something needs to be done about the boundary data!!!
    // For now, we'll assume F satisfies the same boundary conditions that Q
    // does (but this is not true!!)
void SetBndValues(dTensorBC2&, dTensorBC2&);
//SetBndValues(aux, F );

}
