#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "IniParams.h"

using namespace std;

// ------------------------------------------------------------
// Function definitions
void ConSoln( 
    const dTensorBC2& aux,
    const dTensorBC2& q, 
    double t );

// RK functions
void BeforeStep(double dt, dTensorBC2& aux, dTensorBC2& q);
void ConstructL( 
        dTensorBC2& aux,
        dTensorBC2& q,      // setbndy conditions modifies q
        dTensorBC2& Lstar,
        dTensorBC1& smax);
void AfterStep(double dt, dTensorBC2& aux, dTensorBC2& q);

double GetCFL(double dt, double dtmax,
        const dTensorBC2& aux,
        const dTensorBC1& smax);

void BeforeFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);
void AfterFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);
// ------------------------------------------------------------


// ------------------------------------------------------------
// Example function that describes how the main time stepping loop works.
// 
// This routine will terminate the code upon entrance.
// ------------------------------------------------------------
void DogSolveUser(
        dTensorBC2& aux, 
        dTensorBC2& qold,
        dTensorBC2& qnew,
        dTensorBC1& smax,
        double tstart, double tend,int nv, 
        double dtv[], const double cflv[] )
{

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double CFL_max      = cflv[1];  // max   CFL number
    double CFL_target   = cflv[2];  // target CFL number
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const int mx    = qold.getsize(1);
    const int meqn  = qold.getsize(2);
    const int maux  = aux.getsize(2);
    const int mbc   = qnew.getmbc();


    // Allocate storage for this solver
    dTensorBC2   qstar(mx, meqn, mbc);
    dTensorBC2 auxstar(mx, maux, mbc);
    dTensorBC2   Lstar(mx, meqn, mbc);
    dTensorBC2    Lold(mx, meqn, mbc);
    dTensorBC2      q1(mx, meqn, mbc);
    dTensorBC2      q2(mx, meqn, mbc);

    // Set initialize qstar and auxstar values
    qstar.copyfrom( qold );
    auxstar.copyfrom( aux );


    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step = 0;
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            cout << " Error in DogSolveUser.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold
        qold.copyfrom( qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

//          BeforeFullTimeStep(dt, auxstar, aux, qold, qnew);
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

            // do any extra work      
//          AfterFullTimeStep(dt,auxstar,aux,qold,qnew);

            // compute cfl number
            cfl = GetCFL(dt,dtv[2],aux,smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "DogSolve1D ... Step" << setw(5) << n_step;
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
            if (cfl<=CFL_max)
                // accept
            { m_accept = 1; }
            else 
                //reject
            {   
                t = told;
                if( global_ini_params.get_verbosity() )
                {
                    cout<<"FinSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                qold.copyfrom( qnew );
            }

        }

        // compute conservation and print to file
        ConSoln(aux, qnew, t );

    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

}
