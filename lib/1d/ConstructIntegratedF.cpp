#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "assert.h"

// EXPERIMENTAL CODE

void Diff1( double dx, const dTensor2& f, dTensor1& fx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils:
    assert_eq( f.getsize( 2 ), 6 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = (  f.get( m, 1 ) - f.get( m, 6 ) )*(1.0/12.0);
        tmp       += ( -f.get( m, 2 ) + f.get( m, 4 ) )*(2.0/ 3.0);
        fx.set( m, tmp / dx );
    }

}

void Diff2( double dx, const dTensor2& f, dTensor1& fxx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils:
    assert_eq( f.getsize( 2 ), 6 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = ( -f.get( m, 1 ) - f.get( m, 6 ) )*(1.0/12.0);
        tmp       += (  f.get( m, 2 ) + f.get( m, 4 ) )*(4.0/ 3.0);
        tmp       += ( -f.get( m, 3 )                 )*(5.0/ 2.0);
        fxx.set( m, tmp / (dx*dx) );
    }

}

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
void ConstructIntegratedF( double dt, const dTensor2& node, 
    dTensorBC2& aux, dTensorBC2& q,
    dTensorBC1& smax, dTensorBC2& F)
{


    void FluxFunc  (const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void DFluxFunc(const dTensor1& xpts, const dTensor2& Q, const dTensor2& Aux,
	       dTensor3& Dflux);
    void D2FluxFunc(const dTensor1& xpts, const dTensor2& Q, const dTensor2& Aux,
		dTensor4& D2flux);

    // Used for construcing the flux function
    void SampleFunction( 
        int istart, int iend,
        const dTensor2& node,
        const dTensorBC2& qin, 
        const dTensorBC2& auxin,  
        dTensorBC2& Fout,
        void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

    // Problem dimensions (TODO - the boundary data either a) needs one more
    // point, or b) needs to double the number of ghost cells )
    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // Needed to define derivatives
    const double dx    = dogParamsCart1.get_dx();
    const double xlow  = dogParamsCart1.get_xlow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC2 f( mx, meqn, mbc );  // place-holder
    SampleFunction( 1-mbc, mx+mbc, node, q, aux, f, &FluxFunc );

// TODO  - allow for different sized stencils
const int      mpts_sten = 2*mbc;  assert_eq( mpts_sten,      6 );
const int half_mpts_sten =   mbc;  assert_eq( half_mpts_sten, 3 );

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1; i <= mx; i++ )
    {

        // Physical location for this current value:
        dTensor1 xpts( 1 );
        xpts.set( 1, xlow + double(i)*dx - 0.5*dx );

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
        dTensor1 f_t( meqn );
        for( int m1=1; m1 <= meqn; m1++ )
        {
            double tmp = 0.;
            for( int m2=1; m2 <= meqn; m2++ )
            {
                tmp += A.get(1, m1,m2) * fx_val.get(m2);
            }
            f_t.set( m1, tmp );
        }

        // Second-order accuracy:
        for( int m=1; m<=meqn; m++ )
        {
            F.set( i, m, f.get(i,m) - 0.5*dt*f_t.get(m) );
            // F.set( i, m, f.get(i,m) ); // TODO
        }

        // TODO - INCLUDE HIGHER-ORDER TERMS HERE:

    }

    // TODO - something needs to be done about the boundary data!!!
    // For now, we'll assume periodic
void SetBndValues(const dTensor2&, dTensorBC2&, dTensorBC2&);
SetBndValues(node, aux, F );

}



