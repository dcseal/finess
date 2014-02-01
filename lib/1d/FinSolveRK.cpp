#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "RKinfo.h"
#include "FinSolveRK.h"

using namespace std;

void FinSolveRK(
    const dTensor2& node, const dTensor1& prim_vol,      // TODO - remove these params
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir)
{

    // Declare information about the Runge-Kutta method
    const int time_order = dogParams.get_time_order();
    RKinfo rk;
    SetRKinfo(time_order, rk);

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double CFL_max      = cflv[1];  // max   CFL number
    double CFL_target   = cflv[2];  // targe CFL number
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const int mx = qold.getsize(1);
    const int meqn   = qold.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc = qnew.getmbc();

    // Allocate storage for this solver
    dTensorBC2   qstar(mx, meqn, mbc);
    dTensorBC2 auxstar(mx, maux, mbc);
    dTensorBC2   Lstar(mx, meqn, mbc);
    dTensorBC2    Lold(mx, meqn, mbc);
    dTensorBC2      q1(mx, meqn, mbc);
    dTensorBC2      q2(mx, meqn, mbc);

    // Set initialize qstar and auxstar values
    qstar.copyfrom( qold );
    auxstar.copyfrom( aux );

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
            cout << " Error in FinSolveRK.cpp: "<< 
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
            BeforeFullTimeStep(dt, node, prim_vol, auxstar, aux, qold, qnew);

            // Take a full time step of size dt
            switch( time_order )
            {
                case 1:  // First order in time

                    // --------------------------------------------------------
                    // Stage #1 (the only one in this case)
                    rk.mstage = 1;
                    BeforeStep(dt,node,aux,qnew);
                    ConstructL(node, aux, qnew, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,qnew,Lstar,qnew);
                    AfterStep(dt,node,aux,qnew);
                    // --------------------------------------------------------

                    break;

                case 2:  // Second order in time

                    // ---------------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
                    BeforeStep(dt,node,aux,qnew);
                    ConstructL(node,aux,qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,qnew,Lstar,qstar);      
                    AfterStep(dt,node,auxstar,qstar);
                    // ---------------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    BeforeStep(dt,node,auxstar,qstar);
                    ConstructL(node,aux,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,auxstar,qstar,Lstar,qnew);
                    AfterStep(dt,node,aux,qnew); 
                    // ---------------------------------------------------------

                    break;

                case 3:  // Third order in time  (low-storage SSP method)

                    // ---------------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
                    BeforeStep(dt,node,aux,qnew);    
                    ConstructL(node,aux,qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,qnew,Lstar,qstar);
                    AfterStep(dt,node,auxstar,qstar);

                    // ---------------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    BeforeStep(dt,node,auxstar,qstar);
                    ConstructL(node,aux,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,qnew,Lstar,qstar);
                    AfterStep(dt,node,auxstar,qstar);

                    // ---------------------------------------------------------
                    // Stage #3
                    rk.mstage = 3;
                    BeforeStep(dt,node,auxstar,qstar);
                    ConstructL(node,auxstar,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,auxstar,qstar,Lstar,qnew);   
                    AfterStep(dt,node,aux,qnew);
                    // ---------------------------------------------------------

                    break;

                case 4:  // Fourth order in time (10-stages)

                    // -----------------------------------------------
                    //CopyQ(qnew,q1);
                    //CopyQ(q1,q2);
                    q1.copyfrom( qnew );
                    q2.copyfrom( qnew );

                    // Stage: 1,2,3,4, and 5
                    for (int s=1; s<=5; s++)
                    {
                        rk.mstage = s;
                        BeforeStep(dt,node,aux,q1);
                        ConstructL(node,aux,q1,Lstar,smax);
                        if (s==1)
                        {  
                            // CopyQ(Lstar,Lold);  
                            Lold.copyfrom( Lstar );
                        }
                        UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                                rk.beta->get(rk.mstage),dt,node,aux,q1,Lstar,q1);
//                      if (dogParams.using_moment_limiter())
//                      {  ApplyLimiter(node,aux,q1,
//                              &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep(dt,node,aux,q1);
                    }

                    // Temporary storage
                    for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int m=1; m<=meqn; m++)
                    {
                        double tmp = (q2.get(i,m) + 9.0*q1.get(i,m))/25.0;
                        q2.set(i,m, tmp );
                        q1.set(i,m, 15.0*tmp - 5.0*q1.get(i,m) );
                    }

                    // Stage: 6,7,8, and 9
                    for (int s=6; s<=9; s++)
                    {
                        rk.mstage = s;
                        BeforeStep(dt,node,aux,q1);
                        ConstructL(node,aux,q1,Lstar,smax);
                        UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                                rk.beta->get(rk.mstage),dt,node,aux,q1,Lstar,q1);
//                      if (dogParams.using_moment_limiter())
//                      {  ApplyLimiter(node,aux,q1,
//                              &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep(dt,node,aux,q1);
                    }

                    // Stage: 10
                    rk.mstage = 10;
                    BeforeStep(dt,node,aux,q1);
                    ConstructL(node,aux,q1,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,node,aux,q2,Lstar,q1);
//                  if (dogParams.using_moment_limiter())
//                  {  ApplyLimiter(node,aux,q1,
//                          &ProjectRightEig,&ProjectLeftEig);  }
                    AfterStep(dt,node,aux,q1);

                    // CopyQ(q1,qnew);
                    qnew.copyfrom( q1 );
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
                        BeforeStep(dt,node,aux,q1);
                        ConstructL(node,aux,q1,Lstar,smax);
                        if (s==1)
                        {  CopyQ(Lstar,Lold);  }

                        UpdateSoln(
                                rk.gamma->get(1,s), 
                                rk.gamma->get(2,s), 
                                rk.gamma->get(3,s), 
                                rk.delta->get(s), rk.beta->get(s),
                                dt, node, aux, qold, Lstar, q1, q2);

//                      if (dogParams.using_moment_limiter())
//                      {  ApplyLimiter(node,aux,q1,
//                              &ProjectRightEig,&ProjectLeftEig);  }
                        AfterStep(dt,node,aux,q1);
                    }

                    // CopyQ(q1, qnew);
                    qnew.copyfrom( q1 );
                    // -----------------------------------------------          
                    break;

                default:

                    printf("WARNING: torder = %d has not been implemented\n", time_order );

                    break;

            }  // End of switch statement over time-order

            // do any extra work (TODO - add this in later)
            AfterFullTimeStep(dt, node, prim_vol, auxstar, aux, qold, qnew);

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
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln(aux, qnew, t, outputdir);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    DeleteRKInfo(rk);

}
