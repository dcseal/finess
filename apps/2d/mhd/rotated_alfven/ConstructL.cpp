#include <cmath>
#include "assert.h"            // for assert_eq.  Can be removed in future
#include "tensors.h"
#include "dog_math.h"
#include "StateVars.h"
#include "IniParams.h"
#include "ConstructL.h"

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x + g(q,x,t)_y = Psi(q,x,t)
//
void ConstructL( const StateVars& Q, dTensorBC3& Lstar, dTensorBC3& smax)
{

    const dTensorBC3&   q = Q.const_ref_q();
    const dTensorBC3& aux = Q.const_ref_aux();

    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int    mbc = global_ini_params.get_mbc();

    // Size of the WENO stencil
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );

    // The flux, f_{i-1/2, j} and g_{i, j-1/2}.  Recall that the
    // flux lives at the nodal locations, i-1/2, so there is one more term in
    // that vector than on the original grid.  Additionally, for 2D code, there
    // are two terms in the flux: q_t + f_x + g_y = psi.
    dTensorBC3  Fhat(mx+1, my,   meqn, mbc );
    dTensorBC3  Ghat(mx,   my+1, meqn, mbc );

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();

    double alpha1 = 0.;
    double alpha2 = 0.;
    if( global_ini_params.get_global_alpha() )
    {
        // Global wave speed
        GlobalWaveSpd( q, aux, alpha1, alpha2);
    }

    // Normal vector.  This is a carry-over from the DG code.
    dTensor1 nvec(2);

    // --------------------------------------------------------------------- //
    // Compute Fhat{i-1/2, j} - 1st component of the flux function
    // --------------------------------------------------------------------- //
    nvec.set(1, 1.0 );  nvec.set(2, 0.0 );
#pragma omp parallel for
    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,j,m) + q.get(i-1,j,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(maux);
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i-1,j,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1         );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t(         ws+1, maux );

        dTensor2 xvals( ws+1, 2 );
        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int is = i-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( is )*dx - 0.5*dx;
            double yi = ylow + double( j  )*dy - 0.5*dy;
            xvals.set( s, 1, xi );
            xvals.set( s, 2, yi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(is, j, m ) );
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(is, j, ma ) );
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

        // Sample the flux function over the stencil:
        dTensor3 fvals_t( ws+1, meqn, 2 );
        FluxFunc( xvals, qvals_t, auxvals_t, fvals_t );

        // Flux function "f" in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );
        for( int me=1; me <= meqn; me++ )
        for( int s=1; s <= ws+1; s++ )
        {
            f.set(me, s, fvals_t.get( s, me, 1 ) );  // 1st-component - f
            g.set(me, s, fvals_t.get( s, me, 2 ) );  // 2nd-component - g
        }


        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( 1, Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( 1, Auxavg, Qavg,     f, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------

        // -- Compute a local wave speed -- //

        dTensor1 xedge(2), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i-1, j, m) );
            Qr.set(m, q.get(i  , j, m) );
        }

        for( int m=1; m<= maux; m++)
        {
            Auxl.set(m, aux.get(i-1, j, m) );
            Auxr.set(m, aux.get(i  , j, m) );
        }

        // Compute an approximate "fastest" wave speed.
        // TODO - this is redundant in the case of a global value of alpha ...
        // (-DS 6/19/2014)
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( alpha1, Max( abs(s1), abs(s2) ) );
        smax.set( i, j, 1, Max( smax.get(i,j,1), alpha )  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here


        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
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
        ProjectRightEig(1, Auxavg, Qavg, ghat, fhat_loc);
        for( int m=1; m <= meqn; m++ )
        {
            Fhat.set(i, j, m, fhat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Compute Ghat{i, j-1/2} - 2nd-component of the flux function
    // --------------------------------------------------------------------- //
    nvec.set(1, 0.0 );  nvec.set(2, 1.0 );
#pragma omp parallel for
    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,j,m) + q.get(i,j-1,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(maux);
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i,j-1,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1         );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t(         ws+1, maux );
        dTensor2 xvals( ws+1, 2 );
        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int js = j-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( i  )*dx - 0.5*dx;
            double yi = ylow + double( js )*dy - 0.5*dy;
            xvals.set( s, 1, xi );
            xvals.set( s, 2, yi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(i, js, m ) );
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(i, js, ma ) );
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
        dTensor3 fvals_t( ws+1, meqn, 2 );
        FluxFunc( xvals, qvals_t, auxvals_t, fvals_t );

        // Flux function in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );
        for( int me=1; me <= meqn; me++ )
        for( int s=1; s <= ws+1; s++ )
        {
            f.set(me, s, fvals_t.get( s, me, 1 ) );
            g.set(me, s, fvals_t.get( s, me, 2 ) );
        }


        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( 2, Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( 2, Auxavg, Qavg,     g, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------

        // -- Compute a local wave speed -- //

        dTensor1 xedge(2), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i, j-1, m) );
            Qr.set(m, q.get(i, j,   m) );
        }

        for( int m=1; m<= maux; m++)
        {
            Auxl.set(m, aux.get(i, j-1, m) );
            Auxr.set(m, aux.get(i, j,   m) );
        }

        // Compute an approximate "fastest" wave speed.
        // TODO - this is redundant in the case of a global value of alpha ...
        // (-DS 6/19/2014)
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( alpha2, Max( abs(s1), abs(s2) ) );
        smax.set( i, j, 2, Max( smax.get(i,j,2), alpha )  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here

        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
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
        ProjectRightEig(2, Auxavg, Qavg, ghat, fhat_loc);
        for( int m=1; m <= meqn; m++ )
        {
            Ghat.set(i,j,m, fhat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_{i,j} = Lstar = -1/dx( fh_{i+1/2,j} - fh_{i-1/2,j} )
    //                           -1/dy( gh_{i,j+1/2} - gh_{i,j-1/2} ).
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {

        // Compute the source term.
        SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, Lstar, &SourceTermFunc);
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1, j,   m) - Fhat.get(i, j, m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,   j+1, m) - Ghat.get(i, j, m) ) / dy;
                Lstar.set(i,j, m, Lstar.get(i,j,m) + tmp );
            }
        }
    }
    else  // No source term
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1,j,  m) - Fhat.get(i,j,m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,  j+1,m) - Ghat.get(i,j,m) ) / dy;
                Lstar.set(i,j, m, tmp );
            }
        }
    }

    // ---------------------------------------------------------
    // Add extra contributions to Lstar
    // ---------------------------------------------------------
    // LstarExtra(node,aux,q,Lstar);

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
