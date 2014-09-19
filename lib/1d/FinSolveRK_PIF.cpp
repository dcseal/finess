#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "IniParams.h"
#include "RKinfo.h"
#include "FinSolveRK.h"
#include "ConstructL.h"
#include "WenoParams.h"
#include "IniParams.h"

// Used for construcing the flux function
void SampleFunction( 
    int istart, int iend,
    const dTensorBC2& qin, 
    const dTensorBC2& auxin,  
    dTensorBC2& Fout,
    void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

void FluxFunc  (const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&);

double EulerStep( double t, double dt, 
    const dTensorBC2& qold, const dTensorBC2& Lstar,
    dTensorBC2& qnew );

void ConstructL_NOC(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        dTensorBC2& Lstar,
        dTensorBC1& smax);

void ConstructL( 
    const dTensorBC2& aux, 
    const dTensorBC2& q, 
    const dTensorBC2& fstar, 
    dTensorBC2& Lstar, 
    dTensorBC1& smax );

using namespace std;

void DogSolveUser( dTensorBC2& aux, dTensorBC2& qold,
        dTensorBC2& qnew, dTensorBC1& smax,
        double tstart, double tend,int nv, 
        double dtv[], const double cflv[],string outputdir)
{

    void FinSolveRK_PIF(
        dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
        dTensorBC1& smax,
        double tstart, double tend, int nv,
        double dtv[], const double cflv[], string outputdir);
    FinSolveRK_PIF( aux, qold, qnew, smax, tstart, tend, nv,
        dtv, cflv, outputdir);

}

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
void FinSolveRK_PIF(
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir)
{

    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();
    RKinfo rk;
    SetRKinfo(time_order, rk);

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double CFL_max      = cflv[1];  // max   CFL number
    double CFL_target   = cflv[2];  // targe CFL number
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const int mx = qold.getsize(1);
    const int meqn   = qold.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc = qnew.getmbc();


    // Allocate storage for this solver
    dTensorBC2 qstar(mx, meqn, mbc);

    // Flux function (from the PDE)
    dTensorBC2 f1(mx, meqn, mbc);
    dTensorBC2 f2(mx, meqn, mbc);
    dTensorBC2 f3(mx, meqn, mbc);
    dTensorBC2 f4(mx, meqn, mbc);

dTensorBC2 k1(mx, meqn, mbc);
dTensorBC2 k2(mx, meqn, mbc);
dTensorBC2 k3(mx, meqn, mbc);
dTensorBC2 k4(mx, meqn, mbc);

    // Time averaged flux function
    dTensorBC2 fstar(mx, meqn, mbc);

    // Time-averaged right hand side
    dTensorBC2   Lstar(mx, meqn, mbc);

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step = 0;
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step>nv )
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
                    SetBndValues(aux, qnew);
                    SampleFunction( 1-mbc, mx+mbc, qnew, aux, f1, &FluxFunc );

                    // Stage 2:
                    ConstructL_NOC(aux, qnew, Lstar, smax                    );
                    k1.copyfrom( Lstar );
                    t = EulerStep( tn, 0.5*dt, qold, Lstar, qstar            );
                    SetBndValues(aux, qstar);
                    SampleFunction( 1-mbc, mx+mbc, qstar, aux, f2, &FluxFunc );


                    // Stage 3:
                    ConstructL_NOC( aux, qstar, Lstar, smax       );
                    k2.copyfrom( Lstar );
                    t = EulerStep( tn, 0.5*dt, qold, Lstar, qstar );
                    SetBndValues(aux, qstar);
                    SampleFunction( 1-mbc, mx+mbc, qstar, aux, f3, &FluxFunc );


                    // Stage 4:
                    ConstructL_NOC( aux, qstar, Lstar, smax   );
                    k3.copyfrom( Lstar );
                    t = EulerStep( tn, dt, qold, Lstar, qstar );
                    SetBndValues(aux, qstar);
                    SampleFunction( 1-mbc, mx+mbc, qstar, aux, f4, &FluxFunc );

ConstructL_NOC( aux, qstar, Lstar, smax   );
k4.copyfrom( Lstar );

                    // Define right hand side value of q, for final WENO
                    // reconstruction
#pragma omp parallel for
                    for( int i=1-mbc; i <= mx+mbc;   i++ )
                    for( int m=1; m <= meqn; m++ )
                    {
                        double tmp = 0.;
                        tmp +=     f1.get(i,m);
                        tmp += 2.0*f2.get(i,m);
                        tmp += 2.0*f3.get(i,m);
                        tmp +=     f4.get(i,m);
                        fstar.set(i,m, tmp/6.0);
                    }

                    // WENO reconstruction, with projections (for final update)
                    SetBndValues(aux, qnew);
                    ConstructL( aux, qnew, fstar, Lstar, smax );
                    t = EulerStep( tn, dt, qold, Lstar, qnew );
                    SetBndValues(aux, qnew);

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
        SetBndValues(aux, qnew);
        ConSoln(aux, qnew, t, outputdir);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    DeleteRKInfo(rk);

}

// Single Euler Step (on interior points only)
double EulerStep( double t, double dt, 
    const dTensorBC2& qold, const dTensorBC2& Lstar,
    dTensorBC2& qnew )
{

    const int mx      = qnew.getsize(1);     //number of elements in grid
    const int meqn    = qnew.getsize(2);     //number of equations
    const int mbc     = qnew.getmbc();       //number of ghost cells

#pragma omp parallel for
    for( int i=1; i <= mx ; i++ )
    for( int m=1; m <= meqn; m++ )
    {
        double tmp = qold.get(i,m) + dt*Lstar.get(i,m);
        qnew.set(i,m, tmp);
    }

    return t+dt;

}

///////////////////////////////////////////////////////////////////////////////
//
// A ConstructL routine designed intentionally wihtout performing
// characteristic decomposition onto the primitive variables.
//
///////////////////////////////////////////////////////////////////////////////
void ConstructL_NOC(
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
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );

    // The flux, f_{i-1/2}.  Recall that the
    // flux lives at the nodal locations, x_{i-1/2}, so there is one more term in
    // this vector than in q or L.
    dTensorBC2  fhat(mx+1, meqn, mbc );

    // Grid spacing
    const double   xlow = global_ini_params.get_xlow();
    const double     dx = global_ini_params.get_dx();

    double alpha1 = 0.;
    if( global_ini_params.get_global_alpha() )
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

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
//          gp.set( m, s, 0.5*(fvals.get(m,s)      + l_alpha*wvals.get(m,s) )      );
//          gm.set( m, s, 0.5*(fvals.get(m,ws-s+2) - l_alpha*wvals.get(m,ws-s+2) ) );
            gp.set( m, s, 0.5*(fvals.get(m,s)      ) );
            gm.set( m, s, 0.5*(fvals.get(m,ws-s+2) ) );
        }

        // WENO reconstruction
        dTensor2 dGp( meqn, 1 ), dGm( meqn, 1 );
        WenoReconstruct( gp, dGp );
        WenoReconstruct( gm, dGm );

        for( int m=1; m <= meqn; m++ )
        {
            fhat.set(i, m, dGp.get(m,1)+dGm.get(m,1) );
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
    if( global_ini_params.get_source_term() )
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

///////////////////////////////////////////////////////////////////////////////
//
// A routine that performs flux-splitting together with a "time-averaged" flux
// that has already been computed.  That is, fstar is expected to be a
// high-order, in time, approximation to the flux function.
//
///////////////////////////////////////////////////////////////////////////////
void ConstructL( 
    const dTensorBC2& aux, 
    const dTensorBC2& q, 
    const dTensorBC2& fstar, 
    dTensorBC2& Lstar, 
    dTensorBC1& smax )
{

    // Parameters for the current grid (could also use dogParams here)
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    // WENO reconstrution routine
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Size of the WENO stencil
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );

    // The flux, f_{i-1/2}.  Recall that the
    // flux lives at the nodal locations, x_{i-1/2}, so there is one more term in
    // this vector than in q or L.
    dTensorBC2  fhat(mx+1, meqn, mbc );

    // Grid spacing
    const double   xlow = global_ini_params.get_xlow();
    const double     dx = global_ini_params.get_dx();

    double alpha1 = 0.;
    if( global_ini_params.get_global_alpha() )
    {
        // Global wave speed
        GlobalWaveSpd( q, aux, alpha1 );
    }

    // ---------------------------------------------------------------------- //
    // Compute fhat_{i-1/2}
    // ---------------------------------------------------------------------- //
#pragma omp parallel for
    for (int i= 1; i<= mx+1; i++)
    {

        // ------------------------------------------------------------------ //
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ... we should
        //         insert a callback that would permit this to happen.  For
        //         now, we are using simple arithmetic averages.
        // ------------------------------------------------------------------ //
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

        // ------------------------------------------------------------------ //
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // ------------------------------------------------------------------ //

        // Sample q (and f) over the stencil:
        dTensor1 xvals( ws+1 );
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1 );
        dTensor2  fvals( meqn, ws+1 ),  fvals_t ( ws+1, meqn );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t( ws+1, maux );
        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int is = i-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( is )*dx - 0.5*dx;
            xvals.set( s, xi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s,     q.get(is, m ) );
                fvals.set( m, s, fstar.get(is, m ) );
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(is, ma ) );
            }
        }

        // ------------------------------------------------------------------ //
        // The format of Flux and ProjectLeftEig/ProjectRightEig do not
        // contain the same order.  That is, Flux assumes q(1:npts, 1:meqn),
        // whereas the other functions assume q(1:meqn, 1:npts).  For
        // consistency, I will copy back to the latter, because the WENO
        // reconstruction *should* be faster if the list of points is second.
        // (-DS)
        // ------------------------------------------------------------------ //
        ConvertTranspose( qvals,   qvals_t   );
        ConvertTranspose( auxvals, auxvals_t );
        ConvertTranspose( fvals, fvals_t     );

        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( Auxavg, Qavg, fvals, gvals );

        // ------------------------------------------------------------------ //
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // ------------------------------------------------------------------ //

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
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here

        // -- Flux splitting -- //
        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s)      + l_alpha*wvals.get(m,s) )      );
            gm.set( m, s, 0.5*(gvals.get(m,ws-s+2) - l_alpha*wvals.get(m,ws-s+2) ) );
        }

        // ------------------------------------------------------------------ //
        // Part IV: Perform a WENO reconstruction on the characteristic vars.
        // ------------------------------------------------------------------ //
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

    // ---------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_i = Lstar = -1/dx( fh_{i+1/2} - fh_{i-1/2} )
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // ---------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
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

    // ---------------------------------------------------------------------- //
    // Add extra contributions to Lstar
    // ---------------------------------------------------------------------- //
    // LstarExtra(aux,q,Lstar);

}
