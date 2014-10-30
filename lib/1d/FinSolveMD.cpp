#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"
#include "FinSolveMD.h"

void ConstructL( const StateVars& Q, dTensorBC2& Lstar, dTensorBC1& smax);

using namespace std;

// -------------------------------------------------------------------------- //
// Two-stage multiderivative time integration.  Currently, 
// this routine supports a fourth and fifth-order method.
//
// One needs to have a Jacobian defined in DFlux in order to be able to use 
// this routine.
//
// See also: FinSolveRK, FinSolveLxW, FinSolveSDC, and FinSolveUser other solvers.
// -------------------------------------------------------------------------- //
void FinSolveMD( StateVars& Qnew, double tend, double dtv[] )
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

    // Maximum wave speed
    dTensorBC1    smax(mx, mbc);

    // Needed for rejecting a time step
    StateVars Qold( t, mx, meqn, maux, mbc );
    dTensorBC2& qold   = Qold.ref_q();
    dTensorBC2& auxold = Qold.ref_aux();

    // Intermediate storage
    StateVars Qstar( t, mx, meqn, maux, mbc );
    dTensorBC2& qstar   = Qstar.ref_q();
    dTensorBC2& auxstar = Qstar.ref_aux();
    Qstar.copyfrom( Qnew );

    // Right hand side of ODE
    dTensorBC2   Lstar(mx, meqn, mbc);

    // Time-averaged flux function
    dTensorBC2   F(mx, meqn, mbc);

    // Coefficients for the third-order method
    const double A21    = 6.666666666666666e-01;
    const double A31    = 0.0;
    const double A32    = 0.0;

    const double Ahat21 = 2.222222222222222e-01;
    const double Ahat31 = 0.0;
    const double Ahat32 = 0.0;

    const double b1     = 6.250000000000000e-01;
    const double b2     = 3.749999999999999e-01;
    const double b3     = 0.0;

    const double bhat1  = 1.250000000000000e-01;
    const double bhat2  = 1.250000000000000e-01;
    const double bhat3  = 0.0;

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step > nv )
        {
            cout << " Error in DogSolveUser.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold
        Qold.copyfrom( Qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {

            // set current time
            Qnew.set_t( t );
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // ----------------------------------------------------------------
            BeforeFullTimeStep(dt, Qold, Qnew );

            SetBndValues( Qnew  );
            SetBndValues( Qstar );

            switch( global_ini_params.get_time_order() )
            {

// EXPERIMENTAL CASE!
                case 1:

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    A21, Ahat21, aux,     qnew, 
                    0.0, 0.0,  auxstar, qstar,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //

                ConstructIntegratedF( dt, 
                    b1, bhat1, aux, qnew, 
                    b2, bhat2, auxstar, qstar,
                    smax, F);

                // Construct a new right hand side
                ConstructLxWL( aux, qstar, F, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }

                Qnew.set_t( Qnew.get_t() + dt );

                SetBndValues( Qnew  );
                SetBndValues( Qstar );

                // Perform any extra work required:
                AfterStep(dt, Qnew );




break;



                case 3:

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    A21, Ahat21, aux,     qnew, 
                    0.0, 0.0,  auxstar, qstar,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //

                ConstructIntegratedF( dt, 
                    b1, bhat1, aux, qnew, 
                    b2, bhat2, auxstar, qstar,
                    smax, F);

                // Construct a new right hand side
                ConstructLxWL( aux, qstar, F, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }

                Qnew.set_t( Qnew.get_t() + dt );

                SetBndValues( Qnew  );
                SetBndValues( Qstar );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                break;

/*
                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    1.0, 1./3., aux,     qnew, 
                    0.0, 0.0,  auxstar, qstar,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + (2.0*dt/3.0)*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //

//              ConstructIntegratedF( dt, 
//                  4./7., 0., aux, qnew, 
//                  0.0, 0.0, auxstar, qstar,
//                  smax, F);
                ConstructL( Qnew, Lstar, smax);

                // Construct a new right hand side
//              ConstructLxWL( aux, qnew, F, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + (4./7.)*dt*Lstar.vget(k);
                    qnew.vset(k, 7.*tmp/16. );
                }

                ConstructIntegratedF( dt, 
                    0.0, 0.0, aux, qnew, 
                    1., 1./3., auxstar, qstar,
                    smax, F);

                // Construct a new right hand side
                ConstructLxWL( aux, qstar, F, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = ( qstar.vget( k ) + (2.*dt/3.)*Lstar.vget(k) )*9.0/16.0;
                    qnew.vset(k, qnew.vget(k) + tmp );
                }

                Qnew.set_t( Qnew.get_t() + dt );

                SetBndValues( Qnew  );
                SetBndValues( Qstar );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                break;
*/

                case 4:

                // -- Stage 1 -- //
                ConstructIntegratedF( dt, 
                    1.0, 0.25, aux,     qnew, 
                    0.0, 0.0,  auxstar, qstar,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + 0.5*dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    1.0, (1.0/6.0), aux, qnew, 
                    0.0, (1.0/3.0), auxstar, qstar,
                    smax, F);

                // Construct a new right hand side
                ConstructLxWL( aux, qnew, F, Lstar, smax);

                // Update the solution:
  #pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qnew.get_t() + dt );

                SetBndValues( Qnew  );
                SetBndValues( Qstar );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                break;

                case 5:

// Coeffients chosen to optimize region of absolute stability along the
// imaginary axis.
//
// rho = 8.209945182837015e-02 chosen to maximize range of abs. stab. region

                // -- Stage 1 -- //
                ConstructIntegratedF( 2.0/5.0*dt, 
                    1.0, 0.5, 125./8.*8.209945182837015e-02, aux, qnew, 
                    0.0, 0.0, 0.0,                           auxstar, qstar,
                    smax, F);

                ConstructLxWL( aux, qnew, F, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + (2.0/5.0*dt)*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }
                Qstar.set_t( Qnew.get_t() + (2.0/5.0)*dt );

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues( Qnew  );
                SetBndValues( Qstar );

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    1.0, 0.5, (1.0/16.0),     aux, qnew, 
                    0.0, 0.0, (5.0/48.0), auxstar, qstar,
                    smax, F);
                ConstructLxWL( auxstar, qstar, F, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qold.get_t() + dt );

                SetBndValues( Qnew  );
                SetBndValues( Qstar );

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                break;

                default:
                printf("Error.  Time order %d not implemented for multiderivative\n", global_ini_params.get_time_order() );
                exit(1);

            }

            // do any extra work      
            AfterFullTimeStep(dt, Qold, Qnew );

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
                Qnew.copyfrom( Qold  );
            }

        }

        // compute conservation and print to file
        SetBndValues( Qnew );
        ConSoln     ( Qnew );

    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

}
