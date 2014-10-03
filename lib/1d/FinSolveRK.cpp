#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "RKinfo.h"
#include "IniParams.h"
#include "StateVars.h"
#include "FinSolveRK.h"

using namespace std;

void FinSolveRK( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC2& qnew = Qnew.ref_q  ();
    dTensorBC2&  aux = Qnew.ref_aux();

    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();
    RKinfo rk;
    SetRKinfo(time_order, rk);

    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number

    double t            = Qnew.get_t();
    double dt           = dtv[1];   // Start with time step from last frame
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;
    double tmp_t        = 0.;

    const int mx     = qnew.getsize(1);
    const int meqn   = qnew.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = qnew.getmbc();

    // Maximum wave speed
    dTensorBC1    smax(mx, mbc);

    // In case we reject a time step
    StateVars Qold( t, mx, meqn, maux, mbc );
    dTensorBC2& qold   = Qold.ref_q();
    dTensorBC2& auxold = Qold.ref_aux();

    // Intermediate stages
    StateVars Qstar( t, mx, meqn, maux, mbc );
    dTensorBC2&   qstar = Qstar.ref_q();
    dTensorBC2& auxstar = Qstar.ref_aux();

    // Right hand side storage ( for q_t = L( q ) )
    dTensorBC2   Lstar(mx, meqn, mbc);
    dTensorBC2    Lold(mx, meqn, mbc);

    // Local storage (for 4th-order time stepping)
    StateVars    Q1( t, mx, meqn, maux, mbc ); Q1.set_t( Qnew.get_t() );
    dTensorBC2&  q1   = Q1.ref_q();
    dTensorBC2&  aux1 = Q1.ref_aux();

    StateVars    Q2( t, mx, meqn, maux, mbc ); Q2.set_t( Qnew.get_t() );
    dTensorBC2&  q2   = Q2.ref_q();
    dTensorBC2&  aux2 = Q2.ref_aux();

    // Set initialize qstar and auxstar values
    qstar.copyfrom( qold );
    auxstar.copyfrom( aux );

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
            cout << " Error in FinSolveRK.cpp: "<< 
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
            Qnew.set_t( t );
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, Qold, Qnew);

            // Take a full time step of size dt
            switch( time_order )
            {
                case 1:  // First order in time

                    // --------------------------------------------------------
                    // Stage #1 (the only one in this case)
                    rk.mstage = 1;
                    SetBndValues( Qnew  );
                    BeforeStep(dt, Qnew );
                    ConstructL(Qnew, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Qnew);
                    AfterStep(dt, Qnew );
                    // --------------------------------------------------------

                    break;

                case 2:  // Second order in time

                    // ---------------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
                    SetBndValues(  Qnew );
                    BeforeStep(dt, Qnew );
                    ConstructL(Qnew, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Qstar);      
                    AfterStep(dt , Qstar );

                    // ---------------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    SetBndValues(Qstar);
                    BeforeStep(dt, Qstar);
                    ConstructL(Qstar, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qstar, Lstar, Qnew);
                    AfterStep(dt, Qnew );
                    // ---------------------------------------------------------

                    break;

                case 3:  // Third order in time  (low-storage SSP method)

                    // ---------------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
                    SetBndValues(Qnew);
                    BeforeStep(dt,Qnew);    
                    ConstructL(Qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Qstar);
                    AfterStep(dt, Qstar);

                    // ---------------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    SetBndValues(Qstar);
                    BeforeStep(dt, Qstar);
                    ConstructL(Qstar, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Qstar);
                    AfterStep(dt, Qstar);

                    // ---------------------------------------------------------
                    // Stage #3
                    rk.mstage = 3;
                    SetBndValues(Qstar);
                    BeforeStep(dt,Qstar);
                    ConstructL(Qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qstar, Lstar, Qnew);   
                    AfterStep(dt, Qnew);
                    // ---------------------------------------------------------

                    break;

                case 4:  // Fourth order in time (10-stages)

                    // -----------------------------------------------
                    q1.copyfrom( qnew );    Q1.set_t( Qnew.get_t() );
                    q2.copyfrom( qnew );    Q2.set_t( Qnew.get_t() );

                    // Stage: 1,2,3,4, and 5
                    for (int s=1; s<=5; s++)
                    {
                        rk.mstage = s;
                        SetBndValues(Q1);
                        BeforeStep(dt, Q1);
                        ConstructL(Q1, Lstar, smax);
                        if (s==1)
                        {  
                            Lold.copyfrom( Lstar );
                        }
                        UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                                rk.beta->get(rk.mstage), dt, Q1, Lstar, Q1);
                        AfterStep(dt, Q1);
                    }

                    // Temporary storage
                    for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int m=1; m<=meqn; m++)
                    {
                        // TODO - TIME NEEDS TO GET FIXED HERE TOO!
                        double tmp = (q2.get(i,m) + 9.0*q1.get(i,m))/25.0;
                        q2.set(i,m, tmp );
                        q1.set(i,m, 15.0*tmp - 5.0*q1.get(i,m) );
                    }

                    tmp_t = (Q2.get_t() + 9.0*Q1.get_t())/25.0;
                    Q2.set_t( tmp_t );
                    Q1.set_t( 15.0*tmp_t - 5.0*Q1.get_t() );

                    // Stage: 6,7,8, and 9
                    for (int s=6; s<=9; s++)
                    {
                        rk.mstage = s;
                        SetBndValues(Q1);
                        BeforeStep(dt, Q1);
                        ConstructL(Q1, Lstar, smax);
                        UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                                rk.beta->get(rk.mstage), dt, Q1, Lstar, Q1);
                        AfterStep(dt, Q1);
                    }

                    // Stage: 10
                    rk.mstage = 10;
                    SetBndValues(Q1);
                    BeforeStep(dt,Q1);
                    ConstructL(Q1,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Q2, Lstar, Q1);
                    AfterStep(dt, Q1);

                    qnew.copyfrom( q1 );
                    Qnew.set_t( Q1.get_t() );

                    // -----------------------------------------------          
                    break;

                case 5:  // Fifth order in time (8-stages)

                    // -----------------------------------------------
                    // CopyQ(qnew,q1);
                    q1.copyfrom( qnew );
                    q2.setall(0.);

                    for (int s=1; s<=8; s++)
                    {
                        rk.mstage = s;
                        SetBndValues(Q1);
                        BeforeStep(dt,Q1);
                        ConstructL(Q1,Lstar,smax);
                        if (s==1)
                        {  CopyQ(Lstar,Lold);  }

                        UpdateSoln(
                                rk.gamma->get(1,s), 
                                rk.gamma->get(2,s), 
                                rk.gamma->get(3,s), 
                                rk.delta->get(s), rk.beta->get(s),
                                dt,  aux, qold, Lstar, q1, q2);

                        AfterStep(dt, Q1);
                    }
                    qnew.copyfrom( q1 );
                    // -----------------------------------------------          
                    break;

                default:

                    printf("WARNING: torder = %d has not been implemented\n", time_order );

                    break;

            }  // End of switch statement over time-order

            // do any extra work (TODO - add this in later)
            AfterFullTimeStep(dt, Qold, Qnew );

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
                dt    = Min( dtv[2], dt*CFL_target/cfl );
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
                qnew.copyfrom( qold );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        SetBndValues( Qnew );
        ConSoln( Qnew );

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    DeleteRKInfo(rk);

}
