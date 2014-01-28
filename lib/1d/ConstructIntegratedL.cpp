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
    void DFluxFunc (const dTensor1&,const dTensor2&,const dTensor2&,dTensor3&);
    void D2FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor4&);

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
    const int mbc    = q.getmbc();

    // Needed to define derivatives
    const double dx  = dogParamsCart1.get_dx();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC2 f( mx, meqn, mbc );  // place-holder
    SampleFunction( 1-mbc, mx+mbc, node, q, aux, f, &FluxFunc );

// TODO  - allow for different sized stencils
const int      mpts_sten = 2*mbc;  assert_eq( mpts_sten, 6 );
const int half_mpts_sten =   mbc;  assert_eq( mpts_sten, 3 );

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1; i <= mx; i++ )
    {
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
        dTensor2 A( meqn, meqn );
        //DFluxFunc( );

    }

    // TODO - something needs to be done about the boundary data!!!
    // For now, we'll assume periodic
void SetBndValues(const dTensor2&, dTensorBC2&, dTensorBC2&);
SetBndValues(node, aux, F   );
//  SetBndValues(node, aux, f   );
//  SetBndValues(node, aux, fx  );
//  SetBndValues(node, aux, fxx );
//  SetBndValues(node, aux, qx  );

    // Construct the time-integrated flux function at each point:
//  for( int i=1; i <= mx+mbc; i++ )
//  for( int m=1; m <= meqn;   m++ )
//  {
//      double tmp = 

//  }

    // Now that we have the integrated flux function F, we can integrate it:
//  ConstructTI_L( node, aux, q, F, Lstar, smax);

}


// Time-integrated right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
// EXPERIMENTAL CODE
void ConstructTI_L(
        const dTensor2& node,
        const dTensorBC2& aux,
        const dTensorBC2& q,
        const dTensorBC2& F,  // <-- new term: integrated flux
        dTensorBC2& Lstar,
        dTensorBC1& smax)
{

    // Boundary conditions
    //
    // TODO - this should be moved before ConstructL is called (-DS)
//  void SetBndValues(const dTensor2&, dTensorBC2&, dTensorBC2&);
//  SetBndValues(node, aux, q);

    // User supplied functions:
    void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void ProjectLeftEig( const dTensor1& Aux_ave, const dTensor1& Q_ave, 
        const dTensor2& Qvals, dTensor2& Wvals);
    void ProjectRightEig(const dTensor1& Aux_ave, const dTensor1& Q_ave, 
                         const dTensor2& Wvals, dTensor2& Qvals);
    void SetWaveSpd(const dTensor1& xedge, const dTensor1& Ql,
            const dTensor1& Qr, const dTensor1& Auxl, const dTensor1& Auxr,
            double& s1,double& s2);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);

    // Routine for WENO reconstrution
    void WenoReconstruct( const dTensor2& gin, dTensor2& diff_g );

    // Routine to deal with the silly mess where the Fluxes and the
    // Projections are all defined separately.
    void ConvertTranspose( const dTensor2& qin, dTensor2& qout );

    // Parameters for the current grid
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

// TODO - "weno stencil" depends on dogParams.get_space_order(), and ws / 2
// should equal mbc.  This should be added somewhere in the code. 
// (Derived parameters? -DS)
const int  r = 3;  // order = 2*r-1
const int ws = 5;  // Number of points for the weno-reconstruction
assert_eq( mbc, 3 );

    // The flux, f_{i-1/2}.  Recall that the
    // flux lives at the nodal locations, i-1/2, so there is one more term in
    // that vector than on the original grid.
    dTensorBC2  fhat(mx+1, meqn, mbc );

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double   xlow = dogParamsCart1.get_xlow();
    const double     dx = dogParamsCart1.get_dx();

    // ---------------------------------------------------------
    // Compute fhat_{i-1/2}
    // ---------------------------------------------------------
#pragma omp parallel for
    for (int i= 1; i<= mx+1; i++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,m) + q.get(i-1,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(iMax(maux, 1 ) );
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,ma) + aux.get(i-1,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( iMax(maux,1), ws+1         );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t(         ws+1, iMax(maux,1) );
        dTensor1 xvals( ws+1 );
        dTensor2 fvals( meqn, ws+1 );
        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int is = i-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( is )*dx - 0.5*dx;
            xvals.set( s, xi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(is, m ) );
                fvals.set( m, s, F.get(is, m ) ); // <-- NEW part (sample integrated flux)
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(is, ma ) );
            }
        }

        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( Auxavg, Qavg, fvals, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------


        // -- Compute a local wave speed -- //

        dTensor1 xedge(1), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(iMax(1,maux)), Auxr(iMax(1,maux));
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i-1,m) );
            Qr.set(m, q.get(i  ,m) );
        }

        for( int m=1; m<= maux; m++)
        {
            Auxl.set(m, aux.get(i-1,m) );
            Auxr.set(m, aux.get(i  ,m) );
        }

        double s1,s2;
        SetWaveSpd(xedge, Ql, Qr, Auxl, Auxr, s1, s2);  // application specific
        const double alpha = Max( abs(s1), abs(s2) );
        smax.set( i, alpha  );
        const double l_alpha = 1.1*alpha;  // extra safety factor added here

        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s)      + l_alpha*wvals.get(m,s) )      );
            gm.set( m, s, 0.5*(gvals.get(m,ws-s+2) - l_alpha*wvals.get(m,ws-s+2) ) );
        }


        // --------------------------------------------------------------------
        // Part IV: Perform a WENO reconstruction on the characteristic vars.
        // --------------------------------------------------------------------
        dTensor2 dGp( meqn, 1 ), dGm( meqn, 1 );
        WenoReconstruct( gp, dGp );
        WenoReconstruct( gm, dGm );

        // add and convert back to the conserved quantities
        dTensor2 ghat( meqn, 1 );
        for (int m=1; m<=meqn; m++)
        {
            ghat.set(m, 1, dGp.get(m,1) + dGm.get(m,1) );
        }

        dTensor2 fhat_loc( ghat );
        ProjectRightEig(Auxavg, Qavg, ghat, fhat_loc);
        for( int m=1; m <= meqn; m++ )
        {
            fhat.set(i, m, fhat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_i = Lstar = -1/dx( fh_{i+1/2} - fh_{i-1/2} )
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( dogParams.get_source_term() )
    {
        printf("Error: source-term not implemented for Lax-Wendroff method\n");
        exit(1);
        // Compute the source term.
//      SampleFunction( 1-mbc, mx+mbc, node, q, aux, Lstar, &SourceTermFunc);
//  #pragma omp parallel for
//          for (int i=1; i<=mx; i++)
//          {
//              for (int m=1; m<=meqn; m++)
//              {
//                  Lstar.set(i,m, Lstar.get(i,m) -(fhat.get(i+1,m)-fhat.get(i,m))/dx );
//              }
//          }
    }
    else
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        {
            for (int m=1; m<=meqn; m++)
            {
                Lstar.set(i,m, -(fhat.get(i+1,m)-fhat.get(i,m))/dx );
            }
        }
    }
    // ---------------------------------------------------------
    // Add extra contributions to Lstar
    // ---------------------------------------------------------
    // LstarExtra(node,aux,q,Lstar);

}
