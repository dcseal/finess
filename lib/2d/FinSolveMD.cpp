#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "FinSolveLxW.h"     // functions directly called from this function
#include "DogParams.h"
#include "DogParamsCart2.h"

// ------------------------------------------------------------
// Multiderivative integration
//
// These functions are for the two-stage methods.  One contains
// two-derivatives, and the second contains three derivatives.
void ConstructIntegratedR( double dt, 
    double alpha1, double beta1,
    const dTensorBC3& aux1, const dTensorBC3& q1,
    double alpha2, double beta2,
    const dTensorBC3& aux2, const dTensorBC3& q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1, double charlie1,
    const dTensorBC3& aux1, const dTensorBC3& q1,
    double alpha2, double beta2, double charlie2,
    const dTensorBC3& aux2, const dTensorBC3& q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);
// ------------------------------------------------------------


using namespace std;

void FinSolveMD(
    dTensorBC3& aux, dTensorBC3& qold, dTensorBC3& qnew, 
    dTensorBC3& smax,
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

    const double xlow = dogParamsCart2.get_xlow();
    const double ylow = dogParamsCart2.get_ylow();
    const int     mbc = dogParamsCart2.get_mbc();

    const int mx   = dogParamsCart2.get_mx();
    const int my   = dogParamsCart2.get_my();

    const int meqn   = dogParams.get_meqn();
    const int maux   = dogParams.get_maux();

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
    qstar.copyfrom( qnew );   

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
        qstar.copyfrom( qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while( m_accept==0 )
        {

            // set current time
            double told = t;
            dogParams.set_time( t );

            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, aux, aux, qold, qnew);

            SetBndValues(aux,      qnew);
            SetBndValues(auxstar, qstar);

            // ---------------------------------------------------------
            // Take a full time step of size dt
            switch( dogParams.get_time_order() )
            {


                case 4:

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // -- Stage 1 -- //
                ConstructIntegratedR( 0.5*dt, aux, qnew, smax, F, G);

                // That call is equivalent to the following call:
                // Note that the dt has been rescaled in order to retain the
                // correct units for the flux splitting that will occur in a
                // second.
//              ConstructIntegratedR( 0.5*dt, 
//                  1.0, 0.5, aux,     qnew, 
//                  0.0, 0.0, auxstar, qstar,
//                  smax, F, G);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, G, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + 0.5*dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, auxstar, qstar );

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // -- Stage 2 -- //
                ConstructIntegratedR( dt, 
                    1.0, (1.0/6.0), aux, qnew, 
                    0.0, (1.0/3.0), auxstar, qstar,
                    smax, F, G);
                ConstructLxWL( auxstar, qstar, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // Perform any extra work required:
                AfterStep(dt, auxstar, qstar );

                break;

                case 5:

                // Coeffients chosen to optimize region of absolute stability 
                // along the imaginary axis.
                //
                // rho = 8.209945182837015e-02 chosen to maximize range of 
                //                             absolute stability region

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // -- Stage 1 -- //
                ConstructIntegratedR( 2.0/5.0*dt, 
                    1.0, 0.5, 125./8.*8.209945182837015e-02, aux, qnew, 
                    0.0, 0.0, 0.0,                           auxstar, qstar,
                    smax, F, G);

                ConstructLxWL( aux, qnew, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + (2.0/5.0*dt)*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, auxstar, qstar );

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // -- Stage 2 -- //
                ConstructIntegratedR( dt, 
                    1.0, 0.5, (1.0/16.0),     aux, qnew, 
                    0.0, 0.0, (5.0/48.0), auxstar, qstar,
                    smax, F, G);
                ConstructLxWL( auxstar, qstar, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // Perform any extra work required:
                AfterStep(dt, auxstar, qstar );

                break;

                default:
                printf("Error.  Time order %d not implemented for multiderivative\n", dogParams.get_time_order() );
                exit(1);

            }

            // do any extra work
            AfterFullTimeStep(dt, auxstar, aux, qold, qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( dogParams.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveMD2D ... Step" << setw(5) << n_step;
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
                    cout<<"FinSolveMD2D rejecting step...";
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

}
