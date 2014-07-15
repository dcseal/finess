#include <cmath>
#include "assert.h"            // for assert_eq.  Can be removed in future
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "tensors.h"
#include "dog_math.h"
#include "WenoParams.h"
#include "ConstructL.h"

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
// It is expected that the user sets the boundary conditions before calling this
// routine.
void ConstructL(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        dTensorBC2& Lstar,
        dTensorBC1& smax)
{

    // Parameters for the current grid (could also use dogParams here)
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    // WENO reconstrution routine
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Size of the WENO stencil
    const int ws = dogParams.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );

    // The flux, f_{i-1/2}.  Recall that the
    // flux lives at the nodal locations, x_{i-1/2}, so there is one more term in
    // this vector than in q or L.
    dTensorBC2  fhat(mx+1, meqn, mbc );

    // Grid spacing
    const double   xlow = dogParamsCart1.get_xlow();
    const double     dx = dogParamsCart1.get_dx();

    double alpha1 = 0.;
    if( dogParams.get_global_alpha() )
    {
        // Global wave speed
        GlobalWaveSpd( q, aux, alpha1 );
    }

    // ---------------------------------------------------------
    // Compute fhat_{i-1/2}
    // ---------------------------------------------------------
#pragma omp parallel for
    for (int i= 1; i<= mx+1; i++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ... we should
        //         insert a callback that would permit this to happen.  For
        //         now, we are using simple arithmetic averages.
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,m) + q.get(i-1,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(maux);
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,ma) + aux.get(i-1,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1 );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t( ws+1, maux );
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

        // The format of Flux and ProjectLeftEig/ProjectRightEig do not
        // contain the same order.  That is, Flux assumes q(1:npts, 1:meqn),
        // whereas the other functions assume q(1:meqn, 1:npts).  For
        // consistency, I will copy back to the latter, because the WENO
        // reconstruction *should* be faster if the list of points is second.
        // (-DS)
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


        // -- Compute a local wave speed -- //

        dTensor1 xedge(1), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);
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
        const double alpha = Max( alpha1, Max( abs(s1), abs(s2) ) );
        smax.set( i, alpha  );
        const double l_alpha = wenoParams.alpha_scaling*alpha;  // extra safety factor added here

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
        // Compute the source term.
        SampleFunction( 1-mbc, mx+mbc, q, aux, Lstar, &SourceTermFunc);
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        {
            for (int m=1; m<=meqn; m++)
            {
                Lstar.set(i,m, Lstar.get(i,m) -(fhat.get(i+1,m)-fhat.get(i,m))/dx );
            }
        }
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
    // LstarExtra(aux,q,Lstar);

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
