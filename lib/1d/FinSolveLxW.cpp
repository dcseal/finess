#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "FinSolveLxW.h"   // functions directly called from this routine

using namespace std;

void FinSolveLxW(
    const dTensor2& node, const dTensor1& prim_vol,      // TODO - remove these params
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir)
{

    // Declare information about the Runge-Kutta method
    const int time_order = dogParams.get_time_order();

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double CFL_max      = cflv[1];  // max   CFL number
    double CFL_target   = cflv[2];  // targe CFL number
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const int mx     = qold.getsize(1);
    const int meqn   = qold.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = qnew.getmbc();

    // Total number of entries in the vector:
    const int numel = qold.numel();

    // Allocate storage for this solver

    // Flux function, and (linear) derivatives of the flux function
    dTensorBC2    f(mx, meqn, mbc);
    dTensorBC2   fx(mx, meqn, mbc);
    dTensorBC2  fxx(mx, meqn, mbc);

    // Time-integrated flux function:
    dTensorBC2    F(mx, meqn, mbc );

    // Single-derivative of state vector (using conserved variables):
    dTensorBC2   qx(mx, meqn, mbc);

    // Set initialize qstar and auxstar values
    // TODO - we can use the 'copyfrom' routine from the tensor class (-DS)
    //qstar.copyfrom( qold );
    //auxstar.copyfrom( aux );
    dTensorBC2 auxold( mx, meqn, mbc );
    auxold.copyfrom( aux );

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
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, node, prim_vol, aux, aux, qold, qnew);

            // ---------------------------------------------------------
            // Take a full time step of size dt
            BeforeStep(dt,node,aux,qnew);
            ConstructIntegratedL( dt, node, aux, qnew, f, fx, fxx, qx, smax, F);
// Perform a WENO reconstruction on F:  TODO
            
            // Update the solution:
//  #pragma omp parallel for
//          for( int k=0; k < numel; k++ )
//          {
//              double tmp = qnew.vget( k ) + F.vget(k);
//              qnew.vset(k, tmp );
//          }

            // Perform any extra work required:
            AfterStep(dt,node,aux,qnew);
            // ---------------------------------------------------------

            // do any extra work (TODO - add this in later)
            AfterFullTimeStep(dt, node, prim_vol, auxold, aux, qold, qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], prim_vol, aux, smax);

            // output time step information
            if( dogParams.get_verbosity() )
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
                if( dogParams.get_verbosity() )
                {
                    cout<<"FinSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                // CopyQ(qold, qnew);
                qnew.copyfrom( qold );
                auxold.copyfrom( aux );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln(node, aux, qnew, t, outputdir);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

}
