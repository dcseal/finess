#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "RKinfo.h"           // Coefficients for the RK method
#include "FinSolveRK.h"       // Functions directly called from this function
#include "DogParams.h"
#include "DogParamsCart3.h"

using namespace std;

void FinSolveRK(
    dTensorBC4& aux, dTensorBC4& qold, dTensorBC4& qnew, 
    dTensorBC4& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir)
{

    // Declare information about the Runge-Kutta method
    const int time_order = dogParams.get_time_order();
    RKinfo rk;
    SetRKinfo(time_order, rk);

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double CFL_max      = cflv[1];  // max    CFL number
    double CFL_target   = cflv[2];  // target CFL number
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const double xlow = dogParamsCart3.get_xlow();
    const double ylow = dogParamsCart3.get_ylow();
    const double zlow = dogParamsCart3.get_zlow();
    const int     mbc = dogParamsCart3.get_mbc();

    const int mx   = dogParamsCart3.get_mx();
    const int my   = dogParamsCart3.get_my();
    const int mz   = dogParamsCart3.get_mz();

    const int meqn   = dogParams.get_meqn();
    const int maux   = dogParams.get_maux();

    // Allocate storage for this solver
    dTensorBC4   qstar(mx, my, mz, meqn, mbc);
    dTensorBC4 auxstar(mx, my, mz, maux, mbc);
    dTensorBC4   Lstar(mx, my, mz, meqn, mbc);
    dTensorBC4    Lold(mx, my, mz, meqn, mbc);
    dTensorBC4      q1(mx, my, mz, meqn, mbc);
    dTensorBC4      q2(mx, my, mz, meqn, mbc);

    // Set initialize qstar and auxstar values
    // TODO - we can use the 'copyfrom' routine from the tensor class (-DS)
    qstar.copyfrom( qold   );
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
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, auxstar, aux, qold, qnew);

            // Take a full time step of size dt
            switch( time_order )
            {
                case 1:  // First order in time (Forward-Euler)

                    // --------------------------------------------------------
                    // Stage #1 (the only one in this case)
                    rk.mstage = 1;
                    BeforeStep(dt,aux,qnew);
                    ConstructL( aux, qnew, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,aux,qnew,Lstar,qnew);
                    AfterStep(dt,aux,qnew);
                    // --------------------------------------------------------

                    break;

                case 2:  // Second order in time

                    // ---------------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
                    BeforeStep(dt,aux,qnew);
                    ConstructL(aux,qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,aux,qnew,Lstar,qstar);      
                    AfterStep(dt,auxstar,qstar);

                    // ---------------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    BeforeStep(dt,auxstar,qstar);
                    ConstructL(aux,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,auxstar,qstar,Lstar,qnew);
                    AfterStep(dt,aux,qnew); 
                    // ---------------------------------------------------------

                    break;

                case 3:  // Third order in time  (low-storage SSP method)

//     qnew = alpha1 * qstar + alpha2 * qnew + beta * dt * L( qstar )

// alpha1 = 1.0
// alpha2 = 0.0
// beta   = 1.0

                    // ---------------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
dogParams.set_time( told );
                    BeforeStep(dt,aux,qnew);    
                    ConstructL(aux,qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,aux,qnew,Lstar,qstar);
                    AfterStep(dt,auxstar,qstar);

dogParams.set_time( told + dt );
// alpha1 = 0.75
// alpha2 = 0.25
// beta   = 0.25
                    // ---------------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    BeforeStep(dt,auxstar,qstar);
                    ConstructL(aux,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,aux,qnew,Lstar,qstar);
                    AfterStep(dt,auxstar,qstar);

dogParams.set_time( told + (2.0/3.0)*dt );
// alpha1 = 2/3
// alpha2 = 1/3
// beta   = 2/3
                    // ---------------------------------------------------------
                    // Stage #3
                    rk.mstage = 3;
                    BeforeStep(dt,auxstar,qstar);
                    ConstructL(auxstar,qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,auxstar,qstar,Lstar,qnew);   
                    AfterStep(dt,aux,qnew);
                    // ---------------------------------------------------------

                    break;

                case 4:  // Fourth order in time (10-stages)

                    // -----------------------------------------------
                    q1.copyfrom( qnew );
                    q2.copyfrom(   q1 );

                    // Stage: 1,2,3,4, and 5
                    for (int s=1; s<=5; s++)
                    {
                        rk.mstage = s;
                        BeforeStep(dt,aux,q1);
                        ConstructL(aux,q1,Lstar,smax);
                        if (s==1)
                        {  Lold.copyfrom( Lstar ); }
                        UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                                rk.beta->get(rk.mstage),dt,aux,q1,Lstar,q1);
                        AfterStep(dt,aux,q1);
                    }

                    // Temporary storage
                    for(int i = (1-mbc); i <= (mx+mbc); i++)
                    for(int j = (1-mbc); j <= (my+mbc); j++)
                    for(int k = (1-mbc); k <= (mz+mbc); k++)
                    for(int m = 1; m <= meqn; m++)
                    {
                        double tmp = ( q2.get(i,j,k,m) + 9.0*q1.get(i,j,k,m) )/25.0;
                        q2.set(i,j,k,m, tmp );
                        q1.set(i,j,k,m, 15.0*tmp - 5.0*q1.get(i,j,k,m) );
                    }

                    // Stage: 6,7,8, and 9
                    for (int s=6; s<=9; s++)
                    {
                        rk.mstage = s;
                        BeforeStep(dt,aux,q1);
                        ConstructL(aux,q1,Lstar,smax);
                        UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                                rk.beta->get(rk.mstage),dt,aux,q1,Lstar,q1);
                        AfterStep(dt,aux,q1);
                    }

                    // Stage: 10
                    rk.mstage = 10;
                    BeforeStep(dt,aux,q1);
                    ConstructL(aux,q1,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,aux,q2,Lstar,q1);
                    AfterStep(dt,aux,q1);

                    qnew.copyfrom( q1 );
                    // -----------------------------------------------          
                    break;

                case 5:  // Fifth order in time (8-stages)

                    // -----------------------------------------------
                    q1.copyfrom( qnew );
                    q2.setall(0.);

                    for (int s=1; s<=8; s++)
                    {
                        rk.mstage = s;
                        BeforeStep(dt,aux,q1);
                        ConstructL(aux,q1,Lstar,smax);
                        if (s==1)
                        {  Lold.copyfrom(Lstar); }

                        UpdateSoln(
                                rk.gamma->get(1,s), 
                                rk.gamma->get(2,s), 
                                rk.gamma->get(3,s), 
                                rk.delta->get(s), rk.beta->get(s),
                                dt,  aux, qold, Lstar, q1, q2);

                        AfterStep(dt,aux,q1);
                    }

                    qnew.copyfrom( q1 );
                    // -----------------------------------------------          
                    break;

                default:

                    printf("WARNING: torder = %d has not been implemented\n", time_order );

                    break;

            }  // End of switch statement over time-order

            // Do any extra work
            AfterFullTimeStep(dt, auxstar, aux, qold, qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( dogParams.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveRK2D ... Step" << setw(5) << n_step;
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
                    cout<<"FinSolveRK2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                qnew.copyfrom( qold );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln( aux, qnew, t, outputdir);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    DeleteRKInfo(rk);

}
