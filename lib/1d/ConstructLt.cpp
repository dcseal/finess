#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"
#include "StateVars.h"
#include "assert.h"

// Central Finite difference approximations.  See lib/WenoReconstruct.cpp
void Diff1( double dx, const dTensor2& f, dTensor1& fx );
//void Diff2( double dx, const dTensor2& f, dTensor1& fxx );

// !!! NEW TERMS !!! //
void Diff1NC    ( double dx, const dTensor2& f, dTensor1& fx );
//  void Diff2NC    ( double dx, const dTensor2& f, dTensor1& fxx );

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
void ConstructLt( const StateVars& Q, const dTensorBC2& Lstar, dTensorBC1& smax, dTensorBC2& Ltstar)
{

    const dTensorBC2& q   = Q.const_ref_q();
    const dTensorBC2& aux = Q.const_ref_aux();

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    // Placeholder for produce: -A q_t
    dTensorBC2 Aqt( mx, meqn, mbc );

    // TODO  - allow for different sized stencils for different orders (-DS)
    const int mbc_small      = 3;
    const int      mpts_sten = 5;
    const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );

    // Compute -A q_t
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    {

        // Physical location for this current value:
        dTensor1 xpts( 1 );
        xpts.set( 1, xlow + double(i)*dx - 0.5*dx );

        // Compute the Jacobian
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
        DFluxFunc(xpts, q_transpose, a_transpose, A);

        // Compute the product: -A q_t
        dTensor1 f_t( meqn );
        for( int m1=1; m1 <= meqn; m1++ )
        {
            double tmp = 0.;
            for( int m2=1; m2 <= meqn; m2++ )
            {
                tmp += -A.get(1, m1, m2) * Lstar.get(i,m2);
            }
            Aqt.set( i, m1, tmp );
        }

}

    // Compute a (spatial) derivative of -A q_t
    if( global_ini_params.get_space_order() > 1 )
    {
#pragma omp parallel for
    for( int i = 1; i <= mx; i++ )
    {

        // Save the flux function:
        dTensor2 Aqt_vals( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                Aqt_vals.set (m, r, Aqt.get(i+s,m ) );
                s++;
            }
        }

        // Storage for local derivatives:
        dTensor1 Aqt_val( meqn );

        // Compute a FD approximation to the derivatives:
        Diff1( dx, Aqt_vals, Aqt_val );
        //Diff1NC( dx, Aqt_vals, Aqt_val );

        // Second/Third-order accuracy:
        for( int m=1; m<=meqn; m++ )
        {
            Ltstar.set( i, m, Aqt_val.get(m) );
        }

    }
    }
    else
    {
#pragma omp parallel for
    for( int i = 1; i <= mx; i++ )
    {
        // Second-order accuracy
        for( int m=1; m<=meqn; m++ )
        {
            Ltstar.set( i, m, 0.5*(Aqt.get(i+1,m) - Aqt.get(i-1,m))/dx );
//          Ltstar.set( i, m, (Aqt.get(i+1,m) - Aqt.get(i,m))/dx );
//          Ltstar.set( i, m, (Aqt.get(i,m) - Aqt.get(i-1,m))/dx );
        }
                                                          
    }
    }

}
