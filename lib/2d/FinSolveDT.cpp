#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "hooks.h"              // hooks (files that a user may wish to relink)
#include "app_defined.h"        // application (required) files
#include "constructs.h"
#include "misc2d.h"
#include "FinSolveLxW.h"
#include "StateVars.h"
#include "IniParams.h"
using namespace std;


// Stuff used for a single (Euler) step
void SetBndValues( StateVars& Q );
void BeforeStep(double dt, StateVars& Q );
void ConstructL(
        dTensorBC2& aux,
        dTensorBC2& q,      // setbndy conditions modifies q
        dTensorBC2& Lstar,
        dTensorBC1& smax);
void AfterStep(double dt, StateVars& Q );
void ConstructDiffTransformL( const StateVars& Q, const dTensorBC3& F,  dTensorBC3& Lstar, dTensorBC3& smax);
double GetCFL(double dt, double dtmax, const dTensorBC2& aux, const dTensorBC1& smax);
void BeforeFullTimeStep(double dt, const StateVars& Qold, StateVars& Qnew );
void AfterFullTimeStep (double dt, const StateVars& Qold, StateVars& Qnew );

// -------------------------------------------------------------------------- //
//
// FinSolveDT.  Finite difference solver that uses differential transforms.
//
// This solver works with Taylor expansions of the fluxes and computes higher
// derivatives using a recursive definition of the solver.  In order to use
// this solver, the file ConstructDiffTransformL needs to be written locally
// for each application library.  The purpose of putting this file here in the
// main library is to reduce the total amount of code floating around.
// 
//
// See also: FinSolveRK, FinSolveLxW, FinSolveMD, FinSolveUser.  If you don't
// know which solver to use, stick with Runge-Kutta (RK) time stepping.
//
// -------------------------------------------------------------------------- //
void FinSolveUser( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC3& qnew = Qnew.ref_q  ();
    dTensorBC3&  aux = Qnew.ref_aux();

    // Time stepping information
    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number
    double t                  = Qnew.get_t();
    double dt                 = dtv[1];   // Start with time step from last frame
    double cfl                = 0.0;      // current CFL number
    double dtmin              = dt;       // Counters for max and min time step taken
    double dtmax              = dt;

    // Grid information
    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();
    const int meqn = global_ini_params.get_meqn();
    const int maux = global_ini_params.get_maux();
    const int mbc  = global_ini_params.get_mbc();
    const int numel = qnew.numel();

    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

    // Maximum wave speed at each flux interface value.  Note that the size of
    // this tensor is one longer in each direction.
    dTensorBC3 smax( mx+1, my+1, 2, mbc );           

    // Needed for rejecting a time step
    StateVars Qold( t, mx, my, meqn, maux, mbc );
    dTensorBC3& qold   = Qold.ref_q();
    dTensorBC3& auxold = Qold.ref_aux();
    Qold.copyfrom( Qnew );

    // Intermediate stages
    StateVars Qstar( t, mx, my, meqn, maux, mbc );
    dTensorBC3&   qstar = Qstar.ref_q();
    dTensorBC3& auxstar = Qstar.ref_aux();
    Qstar.copyfrom( Qnew );

    // Allocate storage for this solver
    dTensorBC3       F(mx, my, meqn, mbc );  // time-integrated flux
    dTensorBC3       G(mx, my, meqn, mbc );  // time-integrated flux
    dTensorBC3   Lstar(mx, my, meqn, mbc);   // Right hand side of ODE

    // Storage for the MPP limiter
    dTensorBC3* fhat;
    dTensorBC3* fLF;
    dTensorBC3* ghat;
    dTensorBC3* gLF;
    if( global_ini_params.get_mpp_limiter() )
    {
        fhat = new dTensorBC3( mx+1, my, meqn, mbc );
        fLF  = new dTensorBC3( mx+1, my, meqn, mbc );

        ghat = new dTensorBC3( mx, my+1, meqn, mbc );
        gLF  = new dTensorBC3( mx, my+1, meqn, mbc );
    }


    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;                           // Number of time steps taken
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step > nv )
        {
            printf(" Error in FinSolveDT.cpp: "         );
            printf("Exceeded allowed # of time steps \n");
            printf("    n_step = %d\n", n_step          );
            printf("        nv = %d\n", nv              );
            printf("Terminating program.\n"             );
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
            BeforeFullTimeStep(dt, Qold, Qnew);

            // ---------------------------------------------------------
            // Take a full time step of size dt
            BeforeStep(dt, Qnew);
            SetBndValues(Qnew);
            { 
                void ConstructDiffTransformL( double dt, StateVars& Q, dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);
                ConstructDiffTransformL( dt, Qnew, smax, F, G);
                // Construct RHS
                ConstructLxWL( Qnew, F, G, Lstar, smax);

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
            AfterStep(dt, Qnew);
            // ---------------------------------------------------------

            // do any extra work
            AfterFullTimeStep(dt, Qold, Qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveDT2D ... Step" << setw(5) << n_step;
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
                    cout<<"FinSolveDT2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                Qnew.copyfrom( Qold  );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        void ConSoln( const StateVars& Q );
        ConSoln(Qnew);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    // Clean up allocated memory
    if( global_ini_params.get_mpp_limiter() )
    {
        delete fhat;
        delete fLF;
        delete ghat;
        delete gLF;
    }

}
