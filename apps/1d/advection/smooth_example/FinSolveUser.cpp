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

void CheckTotalVariation( const StateVars& Q )
{

    const dTensorBC2& qnew = Q.const_ref_q  ();
    const dTensorBC2&  aux = Q.const_ref_aux();

    // Grid information
    const int mx     = qnew.getsize(1);
    const int meqn   = qnew.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = qnew.getmbc();
    const int numel  = qnew.numel();


    // Total variation
    dTensor1 tv(meqn);
    for (int m=1; m<=meqn; m++)
    {
        tv.set(m, fabs(qnew.get(1,m) - qnew.get(mx,m) ) );
        for (int i=2; i<=mx; i++)
        {
            tv.set(m, tv.get(m) + fabs( qnew.get(i,m) - qnew.get(i-1,m) ) );
        }
    }
    printf("    Total variation = %2.5f\n", tv.get(1) );

}

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
    StateVars Q2( t, mx, meqn, maux, mbc );
    dTensorBC2& q2 = Q2.ref_q();
    dTensorBC2& a2 = Q2.ref_aux();
    Q2.copyfrom( Qnew );

    StateVars Q3( t, mx, meqn, maux, mbc );
    dTensorBC2& q3 = Q3.ref_q();
    dTensorBC2& a3 = Q3.ref_aux();
    Q3.copyfrom( Qnew );

    StateVars Q4( t, mx, meqn, maux, mbc );
    dTensorBC2& q4 = Q4.ref_q();
    dTensorBC2& a4 = Q4.ref_aux();
    Q4.copyfrom( Qnew );

    // Right hand side of ODE
    dTensorBC2   Lstar(mx, meqn, mbc);

    // Time-averaged flux function
    dTensorBC2   F(mx, meqn, mbc);

    // Variables for the Shu-Osher decomposition
    double p21, p31, p32, p41, p42, p43, p51, p52, p53, p54;
    double q21, q31, q32, q41, q42, q43, q51, q52, q53, q54;
    double v[] = {0., 0., 0., 0., 0.};

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

            double r = 0.;
            double rsqd_div_ksqd = 0.;
            switch( global_ini_params.get_time_order() )
            {

                case 3:

                // Shu-Osher coefficients for the third-order method
                p21 = 0.618033988749895;
                p31 = 0.271650617292849;
                p32 = 0.318260723259995;
                q21 = 0.381966011250105;
                q31 = 0.000034591988708;
                q32 = 0.410054067458449;


                // Coefficient from optimal Shu-Osher representation
                r = 1.04;
                rsqd_div_ksqd = r*r / 0.5;

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    p21/r, q21/rsqd_div_ksqd, Qnew, 
                    0.0, 0.0,     Q2,
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
                    p31/r, q31/rsqd_div_ksqd, Qnew, 
                    p32/r, q32/rsqd_div_ksqd, Q2,
                    smax, F);

                // Update the solution:
                ConstructLxWL( a2, q2, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = (p31+q31)*qnew.vget(k) + (p32+q32)*q2.vget(k) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qnew.get_t() + dt );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                SetBndValues(Qnew);
                SetBndValues(Q2);

                break;

                case 4:

                // Shu-Osher coefficients for the fourth-order method
                p21 = 0.618033988749895;
                p31 = 0.362588515112176;
                p32 = 0.207801573327953;
                p41 = 0.144580879241747;
                p42 = 0.110491604448675;
                p43 = 0.426371652664792;

                q21 = 0.381966011250105;
                q31 = 0.;
                q32 = 0.429609911559871;
                q41 = 0.078129569197367;
                q42 = 0.;
                q43 = 0.240426294447419;


                // Coefficient from optimal Shu-Osher representation
                r = 1.392746335264198;
                rsqd_div_ksqd = r*r / 0.5;

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    p21/r, q21/rsqd_div_ksqd, Qnew, 
                    0.0, 0.0, Q2,
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

//              printf("Checking Q2\n");
//              CheckTotalVariation( Q2 );

                SetBndValues(Qnew);
                SetBndValues(Q2);

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    p31/r, q31/rsqd_div_ksqd, Qnew, 
                    p32/r, q32/rsqd_div_ksqd, Q2,
                    smax, F);

                // Update the solution:  (// TODO - use combination of prev. values?)
                ConstructLxWL( a2, q2, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = (p31+q31)*qnew.vget(k) + (p32+q32)*q2.vget(k) + dt*Lstar.vget(k);
                    q3.vset(k, tmp );
                }
                Q3.set_t( Qnew.get_t() + dt );

//              printf("Checking Q3\n");
//              CheckTotalVariation( Q3 );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                SetBndValues(Qnew);
                SetBndValues(Q2);
                SetBndValues(Q3);

                // -- Stage 3 -- //
                ConstructIntegratedF( dt, 
                    p41/r, q41/rsqd_div_ksqd, Qnew, 
                    p42/r, q42/rsqd_div_ksqd, Q2,
                    p43/r, q43/rsqd_div_ksqd, Q3,
                    smax, F);

                // Update the solution:  (// TODO - use combination of prev. values?)
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = (p41+q41)*qnew.vget(k) + 
                        (p42+q42)*q2.vget(k) + (p43+q43)*q3.vget(k) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qnew.get_t() + dt );

//              printf("Checking Qnew\n");
//              CheckTotalVariation( Qnew );

                // Perform any extra work required:
                AfterStep(dt, Qnew );
                SetBndValues(Qnew);

                break;

                case 5:

                // Shu-Osher coefficients for the fifth-order method
                v[0] = 1.0; 
                v[1] = 0.284098093193284;
                v[2] = 0.121597540758194;
                v[3] = 0.000000028591428;
                v[4] = 0.000284064321757;

                // Shu-Osher coefficients for the fourth-order method
                p21 = 0.482803086486156;
                p31 = 0.509529998555284;
                p32 = 0.000000000000004;
                p41 = 0.000000000000000;
                p42 = 0.345844971085767;
                p43 = 0.255534341871175;

                p51 = 0.067790319440499;
                p52 = 0.223408616067964;
                p53 = 0.192058402376403;
                p54 = 0.179899170601193;

                q21 = 0.233098820320560;
                q31 = 0.042689181699210;
                q32 = 0.326183278987308;
                q41 = 0.005301842386952;
                q42 = 0.000003782613078;
                q43 = 0.393315033451601;

                q51 = 0.000000000067575;
                q52 = 0.150752024068163;
                q53 = 0.105203015519228;
                q54 = 0.080604387537218;

                // Coefficient from optimal Shu-Osher representation
                r = 1.354982421031730;
                rsqd_div_ksqd = r*r / 0.5;

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    p21/r, q21/rsqd_div_ksqd, Qnew, 
                    0.0, 0.0, Q2,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    // TODO - fix the notes here (?!)
                    double tmp = (v[1]+p21+q21)*qnew.vget( k ) + dt*Lstar.vget(k);
                    q2.vset(k, tmp );
                }
                Q2.set_t( Qnew.get_t() + dt );

                SetBndValues(Qnew);
                SetBndValues(Q2);

//              printf("Checking Q2\n");
//              CheckTotalVariation( Q2 );

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    p31/r, q31/rsqd_div_ksqd, Qnew, 
                    p32/r, q32/rsqd_div_ksqd, Q2,
                    smax, F);

                // Update the solution:  (// TODO - use combination of prev. values?)
                ConstructLxWL( a2, q2, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = (v[2]+p31+q31)*qnew.vget(k) + (p32+q32)*q2.vget(k) + dt*Lstar.vget(k);
                    q3.vset(k, tmp );
                }
                Q3.set_t( Qnew.get_t() + dt );

//              printf("Checking Q3\n");
//              CheckTotalVariation( Q3 );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                SetBndValues(Qnew);
                SetBndValues(Q2);
                SetBndValues(Q3);

                // -- Stage 3 -- //
                ConstructIntegratedF( dt, 
                    p41/r, q41/rsqd_div_ksqd, Qnew, 
                    p42/r, q42/rsqd_div_ksqd, Q2,
                    p43/r, q43/rsqd_div_ksqd, Q3,
                    smax, F);

                // Update the solution:  (// TODO - use combination of prev. values?)
                // ConstructLxWL( aux, qnew, F, Lstar, smax);
                ConstructLxWL( a3, q3, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = (v[3]+p41+q41)*qnew.vget(k) + 
                        (p42+q42)*q2.vget(k) + (p43+q43)*q3.vget(k) + dt*Lstar.vget(k);
                    q4.vset(k, tmp );
                }
                Q4.set_t( Qnew.get_t() + dt );

                SetBndValues(Qnew);
                SetBndValues(Q2);
                SetBndValues(Q3);
                SetBndValues(Q4);

//              printf("Checking Q4\n");
//              CheckTotalVariation( Q4 );

                // -- Stage 4 -- //
                ConstructIntegratedF( dt, 
                    p51/r, q51/rsqd_div_ksqd, Qnew, 
                    p52/r, q52/rsqd_div_ksqd, Q2,
                    p53/r, q53/rsqd_div_ksqd, Q3,
                    p54/r, q54/rsqd_div_ksqd, Q4,
                    smax, F);

                // Update the solution:  (// TODO - use combination of prev. values?)
                // ConstructLxWL( aux, qnew, F, Lstar, smax);
                ConstructLxWL( a4, q4, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = (v[4]+p51+q51)*qnew.vget(k) + 
                        (p52+q52)*q2.vget(k) + (p53+q53)*q3.vget(k) + (p54+q54)*q4.vget(k) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qnew.get_t() + dt );

//              printf("Checking Qnew\n");
//              CheckTotalVariation( Qnew );

                // Perform any extra work required:
                AfterStep(dt, Qnew );
                SetBndValues(Qnew);

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
