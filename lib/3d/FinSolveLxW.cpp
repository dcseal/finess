#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "FinSolveLxW.h"     // functions directly called from this function
#include "IniParams.h"

using namespace std;

void FinSolveLxW(
    dTensorBC4& aux, dTensorBC4& qold, dTensorBC4& qnew, 
    dTensorBC4& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir)
{

printf("Error: FinSolveLxW not written for 3D code\n");
exit(1);

/*
    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double CFL_max      = cflv[1];  // max   CFL number
    double CFL_target   = cflv[2];  // targe CFL number
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();
    const int     mbc = global_ini_params.get_mbc();

    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();

    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();

    // Total number of entries in the vector:
    const int numel = qold.numel();

    // Allocate storage for this solver

    dTensorBC3 auxstar(mx, my, maux, mbc);   // right hand side (for Euler steps)
    dTensorBC3   qstar(mx, my, meqn, mbc);   // right hand side (for Euler steps)
    dTensorBC3   Lstar(mx, my, meqn, mbc);   // right hand side (for Euler steps)
    dTensorBC3       F(mx, my, meqn, mbc );  // time-integrated flux
    dTensorBC3       G(mx, my, meqn, mbc );  // time-integrated flux

    // Set initialize qstar and auxstar values
    // TODO - we can use the 'copyfrom' routine from the tensor class (-DS)
    if( maux > 0 )
    { auxstar.copyfrom( aux  ); }

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
            cout << " Error in FinSolveLxW.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << "Terminating program." << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold
        // CopyQ(qnew, qold);
        qold.copyfrom( qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while( m_accept==0 )
        {

            // set current time
            double told = t;
            global_ini_params.set_time( t );

            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, aux, aux, qold, qnew);

            // ---------------------------------------------------------
            // Take a full time step of size dt
            BeforeStep(dt, aux, qnew);
            SetBndValues(aux, qnew);
            ConstructIntegratedR( dt, aux, qnew, smax, F, G);

            ConstructLxWL( aux, qnew, F, G, Lstar, smax);  // <-- "new" method

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
            AfterFullTimeStep(dt, auxstar, aux, qold, qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveLxW2D ... Step" << setw(5) << n_step;
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
                    cout<<"FinSolveLxW2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                // CopyQ(qold, qnew);
                qnew.copyfrom( qold  );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln(aux, qnew, t, outputdir);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

*/

}
