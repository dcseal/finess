#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "IniParams.h"
#include "FinSolveLxW.h"   // functions directly called from this routine
#include "StateVars.h"

using namespace std;

// -------------------------------------------------------------------------- //
// Lax-Wendroff time integration.
// this routine supports up to third-order in time.  
//
// One needs to define DFlux and D2FluxFunc in order to be able to use this routine.
//
// See also: FinSolveRK, FinSolveMD, FinSolveSDC, and FinSolveUser other solvers.
// -------------------------------------------------------------------------- //
void FinSolveLxW( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC2& qnew = Qnew.ref_q  ();
    dTensorBC2&  aux = Qnew.ref_aux();

    // Time stepping information
    const double CFL_max        = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target     = global_ini_params.get_desired_cfl();  // target CFL number
    double t                    = Qnew.get_t();
    double dt                   = dtv[1];   // Start with time step from last frame
    double cfl                  = 0.0;      // current CFL number
    double dtmin                = dt;       // Counters for max and min time step taken
    double dtmax                = dt;

    // Grid information
    const int mx     = qnew.getsize(1);
    const int meqn   = qnew.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = qnew.getmbc();
    const int numel  = qnew.numel();
    const double dx  = global_ini_params.get_dx();

    // Maximum wave speed
    dTensorBC1    smax(mx, mbc);

    // Needed for rejecting a time step
    StateVars Qold( t, mx, meqn, maux, mbc );
    dTensorBC2& qold   = Qold.ref_q();
    dTensorBC2& auxold = Qold.ref_aux();

    dTensorBC2  Lstar(mx, meqn, mbc); // right hand side (for Euler steps)
    dTensorBC2    F(mx, meqn, mbc );  // time-integrated flux

    // Storage for the MPP limiter
    dTensorBC2* fhat;
    dTensorBC2* fLF;
    if( global_ini_params.get_mpp_limiter() )
    {
        fhat = new dTensorBC2( mx+1, meqn, mbc );
        fLF  = new dTensorBC2( mx+1, meqn, mbc );
    }

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step>nv )
        {
            cout << " Error in FinSolveLxW.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << "Terminating program." << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold
        Qold.copyfrom( Qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while( m_accept==0 )
        {

            // set current time
            Qnew.set_t( t );
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, Qold, Qnew );

            // ---------------------------------------------------------
            // Take a full time step of size dt
            BeforeStep(dt, Qnew );
            SetBndValues( Qnew  );
            ConstructIntegratedF( dt, aux, qnew, smax, F);

            ConstructLxWL( aux, qnew, F, Lstar, smax); 
            if( global_ini_params.get_mpp_limiter() )
            {

                // Construct the high-order flux
                ConstructLxWL( aux, qnew, F, *fhat, Lstar, smax );
       
                // Construct the low-order flux
                ConstructLFL( Qnew, *fLF );

                // Limit the high-order flux
                ApplyMPPLimiter1D( dt, qnew, *fLF, *fhat );

                // Update the solution:
#pragma omp parallel for
                for( int i=1; i <= mx; i++   )
                for( int m=1; m <= meqn; m++ )
                {
//                  double tmp = (fLF->get(i+1,m)-fLF->get(i,m) );
                    double tmp = (fhat->get(i+1,m)-fhat->get(i,m) );
                    qnew.set(i, m, qnew.get(i,m) - (dt/dx)*tmp );
                }

            }
            else
            { 

                ConstructLxWL( aux, qnew, F, Lstar, smax); 

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }

            }
            Qnew.set_t( Qnew.get_t() + dt );

            // Perform any extra work required:
            AfterStep(dt, Qnew );
            // ---------------------------------------------------------

            // do any extra work
            AfterFullTimeStep(dt, Qold, Qnew );

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
                Qnew.copyfrom( Qold  );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        SetBndValues( Qnew );
        ConSoln     ( Qnew );

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    // Clean up allocated memory
    if( global_ini_params.get_mpp_limiter() )
    {
        delete fhat;
        delete fLF;
    }

}
