#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "IniParams.h"
#include "FinSolveLxW.h"   // functions directly called from this routine

using namespace std;

void FinSolveLxW(
    dTensorBC2& aux, dTensorBC2& qnew, double tstart, 
    double tend, double dtv[] )
{

    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();

    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const int mx     = qnew.getsize(1);
    const int meqn   = qnew.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = qnew.getmbc();

    dTensorBC1    smax(mx, mbc);

    // Total number of entries in the vector:
    const int numel = qnew.numel();

    // Allocate storage for this solver
    dTensorBC2  qold(mx, meqn, mbc);
    dTensorBC2  Lstar(mx, meqn, mbc); // right hand side (for Euler steps)
    dTensorBC2    F(mx, meqn, mbc );  // time-integrated flux

    // Set initialize auxstar values
    dTensorBC2 auxold( mx, maux, mbc );
    if( maux > 0 ){auxold.copyfrom( aux );}

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
        qold.copyfrom( qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while( m_accept==0 )
        {

            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, aux, aux, qold, qnew);

            // ---------------------------------------------------------
            // Take a full time step of size dt
            BeforeStep(dt,aux,qnew);
            SetBndValues(aux, qnew);
            ConstructIntegratedF( dt, aux, qnew, smax, F);
            ConstructLxWL( aux, qnew, F, Lstar, smax);  // <-- "new" method

            // Update the solution:
#pragma omp parallel for
            for( int k=0; k < numel; k++ )
            {
                double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                qnew.vset(k, tmp );
            }

            // Perform any extra work required:
            AfterStep(dt,aux,qnew);
            // ---------------------------------------------------------

            // do any extra work
            AfterFullTimeStep(dt, auxold, aux, qold, qnew);

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
                qnew.copyfrom( qold  );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        SetBndValues(aux, qnew);
        ConSoln(aux, qnew, t );

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

}
