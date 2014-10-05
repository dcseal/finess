#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "RKinfo.h"
#include "FinSolveRK.h"
#include "ConstructL.h"
#include "StateVars.h"
#include "IniParams.h"
#include "app_defined.h"

// Used for construcing the flux function
void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC4& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&));
void FluxFunc(const dTensor2& xpts, const dTensor2& Q,const dTensor2& Aux, dTensor3& flux);

double EulerStep( double t, double dt, 
    const dTensorBC3& qold, const dTensorBC3& Lstar,
    dTensorBC3& qnew );

void ConstructL_NOC( StateVars& Q, dTensorBC3& Lstar, dTensorBC3& smax);

void ConstructL( StateVars& Q,
    dTensorBC3& fstar, dTensorBC3& gstar, 
    dTensorBC3& Lstar, dTensorBC3& smax );

void FinSolveUser( StateVars& Qnew, double tend, double dtv[] )
{
    void FinSolveRK_PIF( StateVars& Qnew, double tend, double dtv[] );
    FinSolveRK_PIF( Qnew, tend, dtv );
}

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// WENO reconstruction without projection onto characteristic variables.  
//
// WARNING: This is an experiemental routine that has not been thoroughly 
//          tested.
//
// Consider stage values, defined by
//
//      k1 = -L( q1 )_{,x}
//      k2 = -L( q2 )_{,x}
//      k3 = -L( q3 )_{,x}
//      k4 = -L( q4 )_{,x}
//
// where L is WENO reconstruction, without projection onto characteristic
// variables.
//
// A bad update would be: q^{n+1} = q^n + dt * ( k1 + 2*(k2+k3) + k4 ).
//
// In this routine, we replace the final derivative with:
//
//      q^{n+1} = q^n - dt * ( f(q1) + 2*( f(q2) + f(q3) ) + f(q4) )_{,x},
//
// where f is the exact flux function, and this derivative is computed using
// WENO, together with the projection onto characteristic variables.
//
///////////////////////////////////////////////////////////////////////////////
void FinSolveRK_PIF( StateVars& Qnew, double tend, double dtv[] )
{

/*
    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();
    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();
    const int meqn = global_ini_params.get_meqn();
    const int maux = global_ini_params.get_maux();
    const int mbc  = global_ini_params.get_mbc();

    dTensorBC3 smax( mx+1, my+1, 2, mbc );           

    // Allocate storage for this solver
    dTensorBC3    qold(mx, my, meqn, mbc);   // Needed for rejecting steps
    dTensorBC3 qstar(mx, my, meqn, mbc);     // Intermediate stage value 

    // Flux function (from the PDE).  These are used to construct a
    // time-averaged flux function.  In the future, a minimal storage method
    // would require performing an ``add-save'' operation on fstar and gstar
    // in place of holding each of these in memory.
    dTensorBC4 R1(mx, my,  meqn, 2, mbc);
    dTensorBC4 R2(mx, my,  meqn, 2, mbc);
    dTensorBC4 R3(mx, my,  meqn, 2, mbc);
    dTensorBC4 R4(mx, my,  meqn, 2, mbc);

    // Right hand side stage values.
    dTensorBC3 k1(mx, my,  meqn, mbc);
    dTensorBC3 k2(mx, my,  meqn, mbc);
    dTensorBC3 k3(mx, my,  meqn, mbc);
    dTensorBC3 k4(mx, my,  meqn, mbc);

    // Time averaged flux function (derived from R1--R4)
    dTensorBC3 fstar(mx, my, meqn, mbc);
    dTensorBC3 gstar(mx, my, meqn, mbc);

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step = 0;
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step > nv )
        {
            cout << " Error in FinSolveRK.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << "Terminating program." << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold
        qold.copyfrom( qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while( m_accept==0 )
        {

            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            double tn = told;
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, aux, aux, qold, qnew);

            // Take a full time step of size dt
            switch( time_order )
            {

                case 4:

                    // Stage 1:
                    SetBndValues(Qnew);
                    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, qnew, aux, R1, &FluxFunc );


                    // Stage 2:
                    ConstructL_NOC(Qnew, k1, smax                    );
                    t = EulerStep( tn, 0.5*dt, qold, k1, qstar            );
                    SetBndValues(Qstar);
                    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, qstar, aux, R2, &FluxFunc );


                    // Stage 3:
                    ConstructL_NOC( Qstar, k2, smax       );
                    t = EulerStep( tn, 0.5*dt, qold, k2, qstar );
                    SetBndValues(Qstar);
                    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, qstar, aux, R3, &FluxFunc );


                    // Stage 4:
                    ConstructL_NOC( Qstar, k3, smax   );
                    t = EulerStep( tn, dt, qold, k3, qstar );
                    SetBndValues(Qstar);
                    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, qstar, aux, R4, &FluxFunc );

                    // Define right hand side value of q, for final WENO
                    // reconstruction
#pragma omp parallel for
                    for( int i=1-mbc; i <= mx+mbc;   i++ )
                    for( int j=1-mbc; j <= my+mbc;   j++ )
                    for( int m=1; m <= meqn; m++ )
                    {
                        double tmp = 0.;
                        tmp +=     R1.get(i,j,m,1);
                        tmp += 2.0*R2.get(i,j,m,1);
                        tmp += 2.0*R3.get(i,j,m,1);
                        tmp +=     R4.get(i,j,m,1);
                        fstar.set(i,j,m, tmp/6.0);

                        tmp = 0.;
                        tmp +=     R1.get(i,j,m,2);
                        tmp += 2.0*R2.get(i,j,m,2);
                        tmp += 2.0*R3.get(i,j,m,2);
                        tmp +=     R4.get(i,j,m,2);
                        gstar.set(i,j,m, tmp/6.0);
                    }

                    // WENO reconstruction, with projections (for final update)
                    SetBndValues(Qnew);
                    ConstructL( Qnew, fstar, gstar, k4, smax );
                    t = EulerStep( tn, dt, qold, k4, qnew );
                    SetBndValues(Qnew);

                    // ------------------------------------------------------------- //
                    // Finished taking a single RK4 time step
                    // ------------------------------------------------------------- //

                    break;

                default:

                    printf("WARNING: torder = %d has not been implemented\n", time_order );
                    break;

            }  // End of switch statement over time-order

            // do any extra work (TODO - add this in later)
            // AfterFullTimeStep(dt, auxstar, aux, qold, qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolve1D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)       // accept
            { m_accept = 1; }
            else                    //reject
            {   
                t = told;
                if( global_ini_params.get_verbosity() )
                {
                    cout<<"FinSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                qnew.copyfrom( qold );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        SetBndValues(Qnew);
        ConSoln(Qnew, t );

    } // End of while loop

    // set initial time step for next call to FinSolve:
    dtv[1] = dt;

*/
}

///////////////////////////////////////////////////////////////////////////////
//
// A routine that performs flux-splitting together with a "time-averaged" flux
// that has already been computed.  That is, fstar and gstar is expected to be a
// high-order, in time, approximation to the flux function.
//
// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x + g(q,x,t)_y = Psi(q,x,t)
//
///////////////////////////////////////////////////////////////////////////////
void ConstructL(
        StateVars& Q,
        dTensorBC3& fstar,      // Time-averaged flux
        dTensorBC3& gstar,      // Time-averaged flux
        dTensorBC3& Lstar,
        dTensorBC3& smax)
{

    dTensorBC3&    q = Q.ref_q  ();
    dTensorBC3&  aux = Q.ref_aux();

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
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1 );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t( ws+1, maux );

        dTensor2  fvals  ( meqn, ws+1  );
        dTensor2  fvals_t( ws+1, meqn  );

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
                qvals.set( m, s, q.get(is, j, m )  );
                fvals.set( m, s, fstar.get(is,j,m) );
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
        ConvertTranspose( fvals,   fvals_t   );
        ConvertTranspose( auxvals, auxvals_t );

//      // Sample the flux function over the stencil:
//      dTensor3 fvals_t( ws+1, meqn, 2 );
//      FluxFunc( xvals, qvals_t, auxvals_t, fvals_t );

        // Flux function "f" in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );
        for( int me=1; me <= meqn; me++ )
        for( int s=1; s <= ws+1; s++ )
        {
            f.set(me, s, fvals_t.get( s, me ) );  // 1st-component - f
//          g.set(me, s, fvals_t.get( s, me, 2 ) );  // 2nd-component - g
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
        dTensor1 Auxl(iMax(1,maux)), Auxr(iMax(1,maux));
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
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1  );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t( ws+1, maux  );

        dTensor2  fvals( meqn, ws+1  ), fvals_t ( ws+1, meqn   );  // <-- NEW

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
                qvals.set( m, s, q.get(i, js,     m ) );
                fvals.set( m, s, gstar.get(i, js, m ) );    // <-- NEW
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
        ConvertTranspose( fvals,   fvals_t   );
        ConvertTranspose( auxvals, auxvals_t );

        // Sample f over the stencil:
//      dTensor3 fvals_t( ws+1, meqn, 2 );
//      FluxFunc( xvals, qvals_t, auxvals_t, fvals_t );

        // Flux function in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );
        for( int me=1; me <= meqn; me++ )
        for( int s=1; s <= ws+1; s++ )
        {
            g.set(me, s, fvals_t.get( s, me ) );
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
        dTensor1 Auxl(iMax(1,maux)), Auxr(iMax(1,maux));
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

///////////////////////////////////////////////////////////////////////////////
//
// A ConstructL routine designed intentionally without performing
// characteristic decomposition onto the primitive variables.
//
///////////////////////////////////////////////////////////////////////////////
void ConstructL_NOC( StateVars& Q,
        dTensorBC3& Lstar, dTensorBC3& smax)
{ 

    dTensorBC3&    q = Q.ref_q  ();
    dTensorBC3&  aux = Q.ref_aux();

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

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1 );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t( ws+1, maux );
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

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
//          gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
//          gm.set( m, s, 0.5*(gvals.get(m,ws-s+2) - l_alpha*wvals.get(m,ws-s+2) ) );
            gp.set( m, s, 0.5*f.get(m,s)      );
            gm.set( m, s, 0.5*f.get(m,ws-s+2) );
        }

        // WENO reconstruction
        dTensor2 dGp( meqn, 1 ), dGm( meqn, 1 );
        WenoReconstruct( gp, dGp );
        WenoReconstruct( gm, dGm );

        for( int m=1; m <= meqn; m++ )
        {
            Fhat.set(i, j, m, dGp.get(m,1) + dGm.get(m,1) );
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

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1 );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t( ws+1, maux );
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

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
//          gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
//          gm.set( m, s, 0.5*(gvals.get(m,ws-s+2) - l_alpha*wvals.get(m,ws-s+2) ) );
            gp.set( m, s, 0.5*(g.get(m,s     ) ) );
            gm.set( m, s, 0.5*(g.get(m,ws-s+2) ) );
        }

        // WENO reconstruction
        dTensor2 dGp( meqn, 1 ), dGm( meqn, 1 );
        WenoReconstruct( gp, dGp );
        WenoReconstruct( gm, dGm );

        for( int m=1; m <= meqn; m++ )
        {
            Ghat.set(i,j,m, dGp.get(m,1)+dGm.get(m,1) );
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

// Single Euler Step (on interior points only)
double EulerStep( double t, double dt, 
    const dTensorBC3& qold, const dTensorBC3& Lstar,
    dTensorBC3& qnew )
{

    const int mx      = qnew.getsize(1);     //number of elements in grid
    const int my      = qnew.getsize(2);     //number of elements in grid
    const int meqn    = qnew.getsize(3);     //number of equations
    const int mbc     = qnew.getmbc();       //number of ghost cells

#pragma omp parallel for
    for( int i=1; i <= mx ; i++ )
    for( int j=1; j <= my ; j++ )
    for( int m=1; m <= meqn; m++ )
    {
        double tmp = qold.get(i,j,m) + dt*Lstar.get(i,j,m);
        qnew.set(i,j,m, tmp);
    }

    return t+dt;

}


