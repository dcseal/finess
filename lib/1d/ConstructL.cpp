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
    void SetWaveSpd(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void SampleFunction( int istart, int iend,
        const dTensor2& node, const dTensorBC2& qin, 
        const dTensorBC2& auxin,  dTensorBC2& Fout,
        void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

    // Routine for WENO reconstrution
    void WenoReconstruct( const dTensor2& gin, dTensor1& diff_g );

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
    dTensor2      qs(mx+1, meqn );
    dTensorBC2    Fp(mx+1, meqn, mbc );
    dTensorBC2    Fm(mx+1, meqn, mbc );
    dTensorBC2  fHat(mx+1, meqn, mbc );
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
    // Part I: Compute Roe Averages
    //         TODO - the User may want to replace this ...
    // ---------------------------------------------------------
//  for( int i=1; i <= mx+1; i++ )
//  for( int m=1; m <= meqn; m++ )
//  {
//      double tmp = 0.5*( q.get(i,m) + q.get(i+1,m) );
//      qs.set(i, m, tmp );
//  }

    // ---------------------------------------------------------
    // Part II: Compute left/right decomposition, f^+ and f^-
    // ---------------------------------------------------------
    for (int i= 1; i<= mx+1; i++)  // TODO - mx+1?
    {

        // Pull the 'left' and 'right' stencils surrounding the point
        // x_{i-1/2}:
        dTensor2 Ql(meqn, ws), Qr(meqn, ws);
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            Ql.set(m, s, q.get( i-1+s-r, m ) );
            Qr.set(m, s, q.get( i+1+r-s, m ) );
        }

        // Convert to characteristic variables
        dTensor2 Wvals(meqn, ws);
        // TODO (this is necessary for a system)
        //      see: ProjectLeftEig.  You'll want to use qs - Q_avg

// Fastest wave speed observed for this element:
//      smax.set(i,   Max(smax_edge, smax.get(i)) );
smax.set( i, 1.0 );  // TODO

dTensor1 dQl( ws );
dTensor1 dQr( ws );
WenoReconstruct( Ql, dQl );
WenoReconstruct( Qr, dQr );

        // Construct fluxes: d/dt q_i = -1/dx( f_{i+1/2} - f_{i-1/2} )
        for (int m=1; m<=meqn; m++)
        {
            Fp.set(i, m, dQl.get(m) );
            Fm.set(i, m, dQr.get(m) );
        }

    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part III: Construct Lstar
    // ---------------------------------------------------------
    for (int i=1; i<=mx; i++)
    {
        for (int m=1; m<=meqn; m++)
        {
            Lstar.set(i,m, -(Fp.get(i+1,m)-Fp.get(i,m))/dx );
        }
    }
    // ---------------------------------------------------------
    // Part IV: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    // LstarExtra(node,aux,q,Lstar);
    // ---------------------------------------------------------

}
