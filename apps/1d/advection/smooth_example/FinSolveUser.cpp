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
// Two-stage multiderivative time integration.
//
// This is an experimental routine that uses a Shu-Osher decomposition of the
// time stepping to update the solution.
//
// See also: FinSolveRK, FinSolveLxW, FinSolveSDC, and FinSolveUser other solvers.
// -------------------------------------------------------------------------- //
void FinSolveUser( StateVars& Qnew, double tend, double dtv[] )
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
    StateVars Q1( t, mx, meqn, maux, mbc );
    dTensorBC2& q1 = Q1.ref_q();
    dTensorBC2& a1 = Q1.ref_aux();
    Q1.copyfrom( Qnew );

    StateVars Q2( t, mx, meqn, maux, mbc );
    dTensorBC2& q2 = Q2.ref_q();
    dTensorBC2& a2 = Q2.ref_aux();
    Q2.copyfrom( Qnew );

    StateVars Q3( t, mx, meqn, maux, mbc );
    dTensorBC2& q3 = Q3.ref_q();
    dTensorBC2& a3 = Q3.ref_aux();
    Q3.copyfrom( Qnew );

    // Right hand side of ODE
    dTensorBC2   Lstar(mx, meqn, mbc);
    dTensorBC2   L1(mx, meqn, mbc);
    dTensorBC2   L2(mx, meqn, mbc);

    // Time-averaged flux function
    dTensorBC2   F(mx, meqn, mbc);
    dTensorBC2  F1(mx, meqn, mbc);
    dTensorBC2  F2(mx, meqn, mbc);

    // Coefficients for the third-order method
    const double p21    = 0.618033988749895;
    const double p31    = 0.271650617292849;
    const double p32    = 0.318260723259995;

    const double q21 = 0.381966011250105;
    const double q31 = 0.000034591988708;
    const double q32 = 0.410054067458449;

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
            SetBndValues( Q1    );
            SetBndValues( Q2    );
            SetBndValues( Q3    );

            double r = 0.;
            double rsqd_div_ksqd = 0.;
            switch( global_ini_params.get_time_order() )
            {

                case 3:

                // Coefficient from optimal Shu-Osher representation
                r = 1.04;
                rsqd_div_ksqd = r*r / 0.5;

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    p21/r, q21/rsqd_div_ksqd, aux, qnew, 
                    0.0, 0.0,     a1, q1,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    q1.vset(k, tmp );
                }
                Q1.set_t( Qnew.get_t() + dt );

                // Perform any extra work required:
                AfterStep(dt, Q1 );

                SetBndValues(Qnew);
                SetBndValues(Q1);

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    p31/r, q31/rsqd_div_ksqd, aux, qnew, 
                    p32/r, q32/rsqd_div_ksqd, a1,    q1,
                    smax, F);

                // Update the solution:
                ConstructLxWL( a1, q1, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = (p31+q31)*qnew.vget(k) + (p32+q32)*q1.vget(k) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qnew.get_t() + dt );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                SetBndValues(Qnew);
                SetBndValues(Q1);

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
