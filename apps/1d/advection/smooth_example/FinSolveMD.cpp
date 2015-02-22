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
// This routine is written as part of the larger scheme to test MSMD methods 
// with and without a Shu-Osher decomposition.  This routine does not use the
// Shu-Osher Decomposition, but FinSolveUser does.
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
    StateVars Q2( t, mx, meqn, maux, mbc );
    dTensorBC2& q2      = Q2.ref_q();
    dTensorBC2& aux2    = Q2.ref_aux();
    Q2.copyfrom( Qnew );

    StateVars Q3( t, mx, meqn, maux, mbc );
    dTensorBC2& q3      = Q3.ref_q();
    dTensorBC2& aux3    = Q3.ref_aux();
    Q3.copyfrom( Qnew );

    // Needed for rejecting a time step
    StateVars Qtmp( t, mx, meqn, maux, mbc );
    dTensorBC2& qtmp   = Qtmp.ref_q();
    dTensorBC2& auxtmp = Qtmp.ref_aux();
    Qtmp.copyfrom( Qnew );

    // Right hand side of ODE
    dTensorBC2   Lstar(mx, meqn, mbc);

    // Time-averaged flux function
    dTensorBC2   F(mx, meqn, mbc);

    // Coefficients for the two-stage, third-order method
//  const double A21    = 0.594243660278112;
//  const double Ahat21 = 0.176562763890364;

//  const double b1     = 0.693990265009012;
//  const double b2     = 0.306009734990988;

//  const double bhat1  = 0.128609262485451; 
//  const double bhat2  = 0.189546392512770;

    // Coefficients for the three-stage, fourth-order method
    const double A21    = 0.443752012194422;
    const double A31    = 0.543193299768317;
    const double A32    = 0.149202742858795;

    const double Ahat21 = 0.098457924163299;
    const double Ahat31 = 0.062758211639901;
    const double Ahat32 = 0.110738910914425;

    const double b1     = 0.515040964378407;
    const double b2     = 0.178821699719783;
    const double b3     = 0.306137335901811;

    const double bhat1  = 0.072864982225864;
    const double bhat2  = 0.073840478463180;
    const double bhat3  = 0.061973770357455;

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
            SetBndValues( Q2    );
            SetBndValues( Q3    );

            switch( global_ini_params.get_time_order() )
            {

                case 3:

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    A21, Ahat21, Qnew, 
                    0.0, 0.0,    Q2,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    q2.vset(k, tmp );
                }
                Q2.set_t( Qnew.get_t() + A21*dt );

                // Perform any extra work required:
                AfterStep(dt, Q2 );

                SetBndValues(Qnew);
                SetBndValues(Q2);

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    b1, bhat1, Qnew, 
                    b2, bhat2, Q2,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qnew.get_t() + dt );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                SetBndValues(Qnew);
                SetBndValues(Q2);

                break;

                case 4:

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    A21, Ahat21, Qnew, 
                    0.0, 0.0,    Q2,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    q2.vset(k, tmp );
                }
                Q2.set_t( Qnew.get_t() + dt );

                // Perform any extra work required:
                AfterStep(dt, Q2 );

                SetBndValues(Qnew);
                SetBndValues(Q2);

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    A31, Ahat31, Qnew, 
                    A32, Ahat32, Q2,
                    smax, F);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = A31*qnew.vget(k) + A32*q2.vget(k);
                    qtmp.vset(k, tmp );
                }
                Qtmp.set_t( Qnew.get_t() + dt );

                // ConstructLxWL( aux, qnew, F, Lstar, smax);
                ConstructLxWL( aux, qtmp, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    q3.vset(k, tmp );
                }
                Q3.set_t( Qnew.get_t() + dt );

                // Perform any extra work required:
                AfterStep(dt, Q3 );

                SetBndValues(Qnew);
                SetBndValues(Q2);
                SetBndValues(Q3);

                // -- Stage 3 -- //
                ConstructIntegratedF( dt, 
                    b1, bhat1, Qnew, 
                    b2, bhat2, Q2,
                    b3, bhat3, Q3,
                    smax, F);

/*
                // Construct a better "approximation" for what to plug into
                // solver
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = b1*qnew.vget(k)+b2*q2.vget(k)+b3*q3.vget(k);
                    qtmp.vset(k, tmp );
                }
                Qtmp.set_t( Qnew.get_t() + dt );


                // Update the solution:
                // ConstructLxWL( aux, qnew, F, Lstar, smax);
                ConstructLxWL( aux, qtmp, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qnew.get_t() + dt );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                SetBndValues(Qnew);
*/
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
