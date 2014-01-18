#include <cmath>
#include "assert.h"            // for assert_eq.  Can be removed in future
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "tensors.h"
#include "dog_math.h"

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
void ConstructL(
        const dTensor2& node,
        dTensorBC2& aux,
        dTensorBC2& q,      // setbndy conditions modifies q
        dTensorBC2& Lstar,
        dTensorBC1& smax)
{

    // Boundary conditions
    //
    // TODO - this should be moved before ConstructL is called (-DS)
    void SetBndValues(const dTensor2&, dTensorBC2&, dTensorBC2&);
    SetBndValues(node, aux, q);

    // User supplied functions:
    void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void ProjectLeftEig( const dTensor1& Aux_ave, const dTensor1& Q_ave, 
        const dTensor2& Qvals, dTensor2& Wvals);
    void ProjectRightEig(const dTensor1& Aux_ave, const dTensor1& Q_ave, 
                         const dTensor2& Wvals, dTensor2& Qvals);
    void SetWaveSpd(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void SampleFunction( int istart, int iend,
        const dTensor2& node, const dTensorBC2& qin, 
        const dTensorBC2& auxin,  dTensorBC2& Fout,
        void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

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
const int ws = 5;
const int  r = 3;  // order = 2*r-1
assert_eq( mbc, 3 );

    // The Roe-average, flux, and source term, respectively.  Recall that the
    // flux lives at the nodal locations, i-1/2, so there is one more term in
    // that vector than on the original grid.
    dTensorBC2  fhat(mx+1, meqn, mbc );
    dTensorBC2   Psi(mx,   meqn, mbc );

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double   xlow = dogParamsCart1.get_xlow();
    const double     dx = dogParamsCart1.get_dx();

    // ---------------------------------------------------------
    // Part 0: compute source term
    // --------------------------------------------------------- 
    Lstar.setall(0.);
//  if( dogParams.get_source_term() )
//  {        
//      // Set source term on computational grid
//      // Set values and apply L2-projection
//      // TODO - replace with SampleFunc
//      SampleFunction( 1-mbc, mx+mbc, node, q, aux, Lstar, &SourceTermFunc);
//  }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part II: Compute left/right decomposition, f^+ and f^-
    // ---------------------------------------------------------
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
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(is, ma ) );
            }
        }
        ConvertTranspose( qvals,   qvals_t   );
        ConvertTranspose( auxvals, auxvals_t );

        // Sample f over the stencil:
        dTensor2 fvals( meqn, ws+1 );  dTensor2 fvals_t( ws+1, meqn );
        FluxFunc( xvals, qvals_t, auxvals_t, fvals_t );
        ConvertTranspose( fvals_t, fvals );

        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( Auxavg, Qavg, fvals, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------
const double alpha = 1.1;  // TODO - insert call to max eigenvalue

/*
printf("printing wvals\n");
for( int s=1; s <= ws+1; s++ )
{
    printf("w = %f\n", wvals.get(1,s) );
}
*/

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s)      + alpha*wvals.get(m,s) )      );
            gm.set( m, s, 0.5*(gvals.get(m,ws-s+2) - alpha*wvals.get(m,ws-s+2) ) );

//printf("    alpha = %f;", alpha );
//printf("    w1, w2 = %2.5e %2.5e; ", wvals.get(m,s), wvals.get(m,ws-s+2) );
//printf("gp, gm = %f, %f\n", gp.get(m,s), gm.get(m,s) );

        }
//printf("\n");
// TODO - Fastest wave speed observed for this element:
smax.set( i, 1.0 );

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

/*
printf("qavg = %2.10e\n", Qavg.get(m) );
printf("ghat = %2.10e\n", ghat.get(m,1) );
printf("fhat = %2.10e\n", fhat.get(i,m) );
*/

        }

    }
    // ---------------------------------------------------------

// Construct fluxes: d/dt q_i = -1/dx( f_{i+1/2} - f_{i-1/2} )

    // ---------------------------------------------------------
    // Part III: Construct Lstar
    // ---------------------------------------------------------
    for (int i=1; i<=mx; i++)
    {
        for (int m=1; m<=meqn; m++)
        {
            Lstar.set(i,m, -(fhat.get(i+1,m)-fhat.get(i,m))/dx );
        }
    }
    // ---------------------------------------------------------
    // Part IV: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    // LstarExtra(node,aux,q,Lstar);
    // ---------------------------------------------------------

}

void ConvertTranspose( const dTensor2& qin, dTensor2& qout )
{
    const int m1 = qin.getsize(1);
    const int m2 = qin.getsize(2);
    assert_eq( m1, qout.getsize(2) );
    assert_eq( m2, qout.getsize(1) );

    for( int i=1; i<= m1; i++ )
    for( int j=1; j<= m2; j++ )
    {
        qout.set(j,i, qin.get(i,j) );
    }

}
