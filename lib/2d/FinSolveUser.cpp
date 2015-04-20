#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "StateVars.h"
#include "IniParams.h"

using namespace std;

// ------------------------------------------------------------
// Function definitions
void ConSoln( const StateVars& Q );

// Stuff used for a single (Euler) step
void SetBndValues( StateVars& Q );
void BeforeStep(double dt, StateVars& Q );
void ConstructL( 
        dTensorBC3& aux,
        dTensorBC3& q,      // setbndy conditions modifies q
        dTensorBC3& Lstar,
        dTensorBC2& smax);
void AfterStep(double dt, StateVars& Q );

double GetCFL(double dt, double dtmax,
        const dTensorBC3& aux,
        const dTensorBC3& smax);

void BeforeFullTimeStep (double dt, StateVars& Qold, StateVars& Q);
void AfterFullTimeStep  (double dt, StateVars& Qold, StateVars& Q);
// ------------------------------------------------------------

// ------------------------------------------------------------
// Example function that describes how the main time stepping loop works.
// 
// This routine will terminate the code upon entrance.
// ------------------------------------------------------------
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

    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();

    // Grid information
    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();
    const int meqn = global_ini_params.get_meqn();
    const int maux = global_ini_params.get_maux();
    const int mbc  = global_ini_params.get_mbc();
    const int numel = qnew.numel();

    // Maximum wave speed at each flux interface value.  Note that the size of
    // this tensor is one longer in each direction.
    dTensorBC3 smax( mx+1, my+1, 2, mbc );           

    // Allocate storage for this solver
    StateVars Qold( t, mx, my, meqn, maux, mbc );
    dTensorBC3& qold   = Qold.ref_q();
    dTensorBC3& auxold = Qold.ref_aux();
    Qold.copyfrom( Qnew );

    // Intermediate stages
    StateVars Qstar( t, mx, my, meqn, maux, mbc );
    dTensorBC3&   qstar = Qstar.ref_q();
    dTensorBC3& auxstar = Qstar.ref_aux();
    Qstar.copyfrom( Qnew );

    dTensorBC3   Lstar(mx, my, meqn, mbc);   // Right hand side of ODE

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;                           // Number of time steps taken
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step       = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step > nv )
        {
            printf(" Error in FinSolveRK.cpp: "         );
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
            // ----------------------------------------------------------------
            //
            //    THIS IS WHERE THE USER-DEFINED TIME-STEPPING SCHEME
            //    SHOULD BE ADDED. IN THE DEFAULT FILE: DogSolveUser.cpp,
            //    THE PROGRAM WILL NOW RETURN AN ERROR MESSAGE.
            // 
            // ----------------------------------------------------------------
            cout << endl;
            cout << " No user-defined time-stepping scheme has been defined yet. " << endl;
            cout << " Copy $FINESS/lib/1d/DogSolveUser.cpp into the current " << endl;
            cout << " directory and modify as needed." << endl << endl;
            exit(1);
            // ----------------------------------------------------------------

            // Do any extra work
            AfterFullTimeStep(dt, Qold, Qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolve2D ... Step" << setw(5) << n_step;
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
                    cout<<"FinSolve2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                Qnew.copyfrom( Qold );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln( Qnew );

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

}
