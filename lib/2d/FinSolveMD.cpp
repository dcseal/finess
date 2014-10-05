#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "FinSolveLxW.h"     // functions directly called from this function
#include "constructs.h" 
#include "StateVars.h"
#include "IniParams.h"

// ------------------------------------------------------------
// Multiderivative integration
//
// These functions are for the two-stage methods.  One contains
// two-derivatives, and the second contains three derivatives.
void ConstructIntegratedR( double dt, 
    double alpha1, double beta1,
    const StateVars& Q1,
    double alpha2, double beta2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1, double charlie1,
    const StateVars& Q1,
    double alpha2, double beta2, double charlie2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);
// ------------------------------------------------------------


using namespace std;

void FinSolveMD( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC3& qnew = Qnew.ref_q  ();
    dTensorBC3&  aux = Qnew.ref_aux();

    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();

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
    dTensorBC3   Lstar(mx, my, meqn, mbc);   // right hand side (for Euler steps)
    dTensorBC3       F(mx, my, meqn, mbc );  // time-integrated flux
    dTensorBC3       G(mx, my, meqn, mbc );  // time-integrated flux

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

            // ---------------------------------------------------------
            // Take a full time step of size dt
            switch( global_ini_params.get_time_order() )
            {


                case 4:

                SetBndValues(Qnew);

                // -- Stage 1 -- //
                ConstructIntegratedR( 0.5*dt, Qnew, smax, F, G);

                // That call is equivalent to the following call:
                // Note that the dt has been rescaled in order to retain the
                // correct units for the flux splitting that will occur in a
                // second.
//              ConstructIntegratedR( 0.5*dt, 
//                  1.0, 0.5, aux,     qnew, 
//                  0.0, 0.0, auxstar, qstar,
//                  smax, F, G);

                // Update the solution:
                ConstructLxWL( Qnew, F, G, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + 0.5*dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }
                Qstar.set_t( Qnew.get_t() + 0.5*dt );

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //
                ConstructIntegratedR( dt, 
                    1.0, (1.0/6.0), Qnew, 
                    0.0, (1.0/3.0), Qstar,
                    smax, F, G);
                ConstructLxWL( Qstar, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qold.get_t() + dt );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                break;

                case 5:

                // Coeffients chosen to optimize region of absolute stability 
                // along the imaginary axis.
                //
                // rho = 8.209945182837015e-02 chosen to maximize range of 
                //                             absolute stability region

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 1 -- //
                ConstructIntegratedR( 2.0/5.0*dt, 
                    1.0, 0.5, 125./8.*8.209945182837015e-02, Qnew, 
                    0.0, 0.0, 0.0,                           Qstar,
                    smax, F, G);

                ConstructLxWL( Qnew, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + (2.0/5.0*dt)*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }
                Qstar.set_t( Qnew.get_t() + (2.0/5.0*dt) );

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //
                ConstructIntegratedR( dt, 
                    1.0, 0.5, (1.0/16.0), Qnew, 
                    0.0, 0.0, (5.0/48.0), Qstar,
                    smax, F, G);
                ConstructLxWL( Qstar, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qold.get_t() + dt );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                break;

                default:
                printf("Error.  Time order %d not implemented for multiderivative\n", global_ini_params.get_time_order() );
                exit(1);

            }

            // do any extra work
            AfterFullTimeStep(dt, Qold, Qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
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
                if( global_ini_params.get_verbosity() )
                {
                    cout<<"FinSolveMD2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                Qnew.copyfrom( Qold  );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln(Qnew);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

}
