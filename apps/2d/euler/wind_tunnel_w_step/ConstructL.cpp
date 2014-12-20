#include <cmath>
#include "assert.h"            // for assert_eq.  Can be removed in future
#include "IniParams.h"
#include "tensors.h"
#include "dog_math.h"
#include "StateVars.h"

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x + g(q,x,t)_y = Psi(q,x,t)
//
void ConstructL(
        StateVars& Q,      // SetBndValues modifies q and aux
        dTensorBC3& Lstar,
        dTensorBC3& smax)
{

    dTensorBC3& q   = Q.ref_q();
    dTensorBC3& aux = Q.ref_aux();

    // Boundary conditions
    //
//  void SetBndValues (dTensorBC3& aux, dTensorBC3& q);  // The "wrong" one.
    void SetBndValuesX(StateVars& Q);  // Only set conditions along x-direction
    void SetBndValuesY(StateVars& Q);  // Only set conditions along y-direction

    // --- User supplied functions --- //
    void FluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor3& flux);
    void ProjectLeftEig( int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
        const dTensor2& Qvals, dTensor2& Wvals);
    void ProjectRightEig(int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
                         const dTensor2& Wvals, dTensor2& Qvals);
    void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge, 
        const dTensor1& Ql,   const dTensor1& Qr, 
        const dTensor1& Auxl, const dTensor1& Auxr,
        double& s1,double& s2);
    void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, 
                const dTensor2& auxvals, dTensor2& source);

    void SampleFunction( 
        int istart, int iend,
        int jstart, int jend,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  dTensorBC3& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Routine to deal with the silly mess where the Fluxes and the
    // Projections are all defined separately.
    void ConvertTranspose( const dTensor2& qin, dTensor2& qout );

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


    // --------------------------------------------------------------------- //
    // Compute Fhat{i-1/2, j} - 1st component of the flux function
    // --------------------------------------------------------------------- //
    SetBndValuesX(Q);

    // Normal vector.  This is a carry-over from the DG code.
    dTensor1 nvec(2);
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
//      for( int m=1; m <= meqn; m++ )
//      {
//          double tmp = 0.5*( q.get(i,j,m) + q.get(i-1,j,m) );
//          Qavg.set(m, tmp );
//      }

        // Densities
        double rhol     = fabs( q.get(i-1,j,1) );
        double sq_rhol  = sqrt( rhol );

        double rhor     = fabs( q.get(i  ,j,1) );
        double sq_rhor  = sqrt( rhor );

        double rho_tmp  = sq_rhor + sq_rhol + 1e-13;


        // Velocities
        double u1l = q.get(i-1,j,2)/rhol;
        double u1r = q.get(i  ,j,2)/rhor;

        double u2l = q.get(i-1,j,3)/rhol;
        double u2r = q.get(i  ,j,3)/rhor;

        double u1m = ( u1l*sq_rhol + u1r*sq_rhor ) / rho_tmp;
        double u2m = ( u2l*sq_rhol + u2r*sq_rhor ) / rho_tmp;

        // Energy and Pressure
        double energyl = q.get(i-1,j,5);
        double energyr = q.get(i,j,5);

        const double gamma = global_ini_params.get_gamma();
        const double gm1   = gamma-1.0;

        double pressl  = gm1*(energyl-0.5e0*rhol*(u1l*u1l+u2l*u2l));
        double pressr  = gm1*(energyr-0.5e0*rhor*(u1r*u1r+u2r*u2r));

        // Enthalpy
        double Hl = (energyl + pressl)/rhol;
        double Hr = (energyr + pressr)/rhor;
        double Hm = ( sq_rhol*Hl + sq_rhor*Hr ) / rho_tmp;
        
        // sound speed
        double cl = sqrt( gamma*pressl / rhol );
        double cr = sqrt( gamma*pressr / rhor );
        double cm = ( gm1*(Hm - 0.5*( u1m*u1m + u2m*u2m ) ) );

        // Now, convert back to "conserved" variable
        double rhom     = sq_rhol*sq_rhor;
        double pressm   = rhom*cm*cm/gamma;
        double Em       = rhom*Hm - pressm;

        Qavg.set( 1, rhom       );
        Qavg.set( 2, rhom * u1m );
        Qavg.set( 3, rhom * u2m );
        Qavg.set( 4, 0.0        );
        Qavg.set( 5, Em         );

        // --------------------------------------------------------------------
        // End of computing Roe averages
        // --------------------------------------------------------------------

        dTensor1 Auxavg(iMax(maux, 1 ) );
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i-1,j,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( iMax(maux,1), ws+1         );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t(         ws+1, iMax(maux,1) );

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

        dTensor1 xedge(1), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(iMax(1,maux)), Auxr(iMax(1,maux));
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );

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
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( abs(s1), abs(s2) );
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
    SetBndValuesY(Q);
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
//      for( int m=1; m <= meqn; m++ )
//      {
//          double tmp = 0.5*( q.get(i,j,m) + q.get(i,j-1,m) );
//          Qavg.set(m, tmp );
//      }

        // Densities
        double rhol     = fabs( q.get(i-1,j,1) );
        double sq_rhol  = sqrt( rhol );

        double rhor     = fabs( q.get(i  ,j,1) );
        double sq_rhor  = sqrt( rhor );

        double rho_tmp  = sq_rhor + sq_rhol + 1e-13;

        // Velocities
        double u1l = q.get(i-1,j,2)/rhol;
        double u1r = q.get(i  ,j,2)/rhor;

        double u2l = q.get(i-1,j,3)/rhol;
        double u2r = q.get(i  ,j,3)/rhor;

        double u1m = ( u1l*sq_rhol + u1r*sq_rhor ) / rho_tmp;
        double u2m = ( u2l*sq_rhol + u2r*sq_rhor ) / rho_tmp;

        // Energy and Pressure
        double energyl = q.get(i-1,j,5);
        double energyr = q.get(i,j,5);

        const double gamma = global_ini_params.get_gamma();
        const double gm1   = gamma-1.0;

        double pressl  = gm1*(energyl-0.5e0*rhol*(u1l*u1l+u2l*u2l));
        double pressr  = gm1*(energyr-0.5e0*rhor*(u1r*u1r+u2r*u2r));

        // Enthalpy
        double Hl = (energyl + pressl)/rhol;
        double Hr = (energyr + pressr)/rhor;
        double Hm = ( sq_rhol*Hl + sq_rhor*Hr ) / rho_tmp;
        
        // sound speed
        double cl = sqrt( gamma*pressl / rhol );
        double cr = sqrt( gamma*pressr / rhor );
        double cm = ( gm1*(Hm - 0.5*( u1m*u1m + u2m*u2m ) ) );

        // Now, convert back to "conserved" variable
        double rhom     = sq_rhol*sq_rhor;
        double pressm   = rhom*cm*cm/gamma;
        double Em       = rhom*Hm - pressm;

        Qavg.set( 1, rhom       );
        Qavg.set( 2, rhom * u1m );
        Qavg.set( 3, rhom * u2m );
        Qavg.set( 4, 0.0        );
        Qavg.set( 5, Em         );

        // --------------------------------------------------------------------
        // End of computing Roe averages
        // --------------------------------------------------------------------

        dTensor1 Auxavg(iMax(maux, 1 ) );
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i,j-1,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( iMax(maux,1), ws+1         );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t(         ws+1, iMax(maux,1) );
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

        dTensor1 xedge(1), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(iMax(1,maux)), Auxr(iMax(1,maux));
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
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
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( abs(s1), abs(s2) );
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

    // Reset smax so it is zero inside the wedge:
    const int istep = (mx/5);
    const int jstep = (my/5)+1;
    for (int i=istep+1; i<= mx+mbc;  i++)
    for (int j=1-mbc;   j<= jstep-1; j++)
    {
        smax.set(i, j, 1, 0.);
        smax.set(i, j, 2, 0.);
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
