#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "RKinfo.h"             // Coefficients for the RK method
#include "FinSolveRK.h"         // Functions directly called from this function
#include "app_defined.h"
#include "StateVars.h"          // Information for state variables
#include "IniParams.h"          // Global parameters accessor

#include "ConstructHJ_L.h"
#include "app_util.h"

using namespace std;

void FinSolveRK( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC3& qnew = Qnew.ref_q  ();
    dTensorBC3&  aux = Qnew.ref_aux();

    // Time stepping information
    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number
    double t                  = Qnew.get_t();
    double dt                 = dtv[1];   // Start with time step from last frame
    double cfl                = 0.0;      // current CFL number
    double dtmin              = dt;       // Counters for max and min time step taken
    double dtmax              = dt;

    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();
    RKinfo rk;
    SetRKinfo(time_order, rk);
    double tmp_t = 0.;                    // used for fourth-order stepping

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

    dTensorBC3   Lstar(mx, my, meqn, mbc);   // Right hand side of ODE
    dTensorBC3    Lold(mx, my, meqn, mbc);
    
    dTensorBC3 Lauxstar(mx, my, maux, mbc);

    // Local storage (for 4th- and 5th-order time stepping)
    StateVars    Q1( t, mx, my, meqn, maux, mbc );
    dTensorBC3&  q1   = Q1.ref_q();
    dTensorBC3&  aux1 = Q1.ref_aux();
    Q1.copyfrom( Qnew );

    StateVars    Q2( t, mx, my, meqn, maux, mbc );
    dTensorBC3&  q2   = Q2.ref_q();
    dTensorBC3&  aux2 = Q2.ref_aux();
    Q2.copyfrom( Qnew );

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;                           // Number of time steps taken
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step       = n_step + 1;

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

//     qnew = alpha1 * qstar + alpha2 * qnew + beta * dt * L( qstar )


                    // ---------------------------------------------------------
                    // Stage #1
                    //      alpha1 = 1.0
                    //      alpha2 = 0.0
                    //      beta   = 1.0
                    // ---------------------------------------------------------
                    rk.mstage = 1;
                    SetBndValues(Qnew);
                    BeforeStep(dt,Qnew);    
                    ConstructL(Qnew,Lstar,smax);
                    ConstructHJ_L(Qnew, Lauxstar);
                    UpdateSolnWithAux(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Lauxstar, Qstar);
                    AfterStep(dt, Qstar);


                    // ---------------------------------------------------------
                    // Stage #2
                    //      alpha1 = 0.75
                    //      alpha2 = 0.25
                    //      beta   = 0.25
                    // ---------------------------------------------------------

                    rk.mstage = 2;
                    SetBndValues(Qstar);
                    BeforeStep(dt, Qstar);
                    ConstructL(Qstar, Lstar, smax);
                    ConstructHJ_L(Qstar, Lauxstar);
                    UpdateSolnWithAux(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Lauxstar, Qstar);
                    AfterStep(dt, Qstar);

                    // ---------------------------------------------------------
                    // Stage #3
                    //      alpha1 = 2/3
                    //      alpha2 = 1/3
                    //      beta   = 2/3
                    // ---------------------------------------------------------

                    rk.mstage = 3;
                    SetBndValues(Qstar);
                    BeforeStep(dt,Qstar);
                    ConstructL(Qstar,Lstar,smax);
                    ConstructHJ_L(Qstar, Lauxstar);
                    UpdateSolnWithAux(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qstar, Lstar, Lauxstar, Qnew);   
                    AfterStep(dt, Qnew);
                    // ---------------------------------------------------------

                    break;

                case 4: // Fourth order in time (10-stages) See Pseudocode 3 in
                        //
                        // "Highly Efficient Strong Stability Preserving Runge-Kutta Methods with
                        // Low-Storage Implementations," David I. Ketcheson, SIAM Journal on Scientific 
                        // Computing, 30(4):2113-2136 (2008)
                        //

                    // -----------------------------------------------
                    Q1.copyfrom( Qnew );
                    Q2.copyfrom( Qnew );

                    // Stage: 1,2,3,4, and 5
                    for (int s=1; s<=5; s++)
                    {
                        rk.mstage = s;
                        SetBndValues(Q1);
                        BeforeStep(dt, Q1);
                        ConstructL(Q1, Lstar, smax);
                        if (s==1)
                        {  Lold.copyfrom( Lstar ); }
                        UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                                rk.beta->get(rk.mstage), dt, Q1, Lstar, Q1);
                        AfterStep(dt, Q1);
                    }

                    // Temporary storage
                    #pragma omp parallel for
                    for( int k=0; k < numel; k++ )
                    {
                        double tmp = (q2.vget(k) + 9.0*q1.vget(k))/25.0;
                        q2.vset(k, tmp );
                        q1.vset(k, 15.0*tmp - 5.0*q1.vget(k) );

                    }

                    if( maux > 0 )
                    {
                        int numel_aux = aux1.numel();
                        #pragma omp parallel for
                        for( int k=0; k < numel_aux; k++ )
                        {
                            double tmp = (aux2.vget(k) + 9.0*aux1.vget(k))/25.0;
                            aux2.vset(k, tmp );
                            aux1.vset(k, 15.0*tmp - 5.0*aux1.vget(k) );

                        }
                    }

                    // Swap the time values as well (L=1)
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

                    Qnew.copyfrom( Q1 );

                    break;

                case 5: // Fifth order in time (8-stages)
                        // TODO - what paper did these coefficients come from?

//                  Q1.copyfrom( Qnew );   // we can remove two replacements
                                           // here
                    q2.setall(0.);      
                    Q2.set_t( 0.);

                    for (int s=1; s<=8; s++)
                    {
                        rk.mstage = s;
                        SetBndValues( Qnew );
                        BeforeStep(dt, Qnew );
                        ConstructL(Qnew, Lstar, smax);
                        if( s==1 )
                        {  Lold.copyfrom( Lstar ); }

                        UpdateSoln( rk.gamma->get(1,s), rk.gamma->get(2,s), rk.gamma->get(3,s), 
                                rk.delta->get(s), rk.beta->get(s), dt,  Qold, Lstar, Qnew, Q2);

                        AfterStep(dt, Qnew);
                    }

// TODO - the time information for this isn't working correctly.
//printf("Q, t1, t2 = %f, %f, %f \n", qnew.get(1,1), Qnew.get_t(), Q2.get_t() );
//assert_lt( fabs( Qnew.get_t() - t ), 1e-8 );

                    break;

                default:

                    printf("WARNING: torder = %d has not been implemented\n", time_order );
                    break;

            }  // End of switch statement over time-order

            // Do any extra work
            AfterFullTimeStep(dt, Qold, Qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
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
                if( global_ini_params.get_verbosity() )
                {
                    cout<<"FinSolveRK2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                Qnew.copyfrom( Qold );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln( Qnew );

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    DeleteRKInfo(rk);

}
