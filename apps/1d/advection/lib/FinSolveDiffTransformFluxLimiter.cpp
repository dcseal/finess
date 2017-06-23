#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"
#include "FinSolveLxW.h"

using namespace std;

// ------------------------------------------------------------
// Function definitions
void ConSoln( const StateVars& Q );

// Stuff used for a single (Euler) step
void SetBndValues( StateVars& Q );
void BeforeStep(double dt, StateVars& Q );
void ConstructL( 
        dTensorBC2& aux,
        dTensorBC2& q,      // setbndy conditions modifies q
        dTensorBC2& Lstar,
        dTensorBC1& smax);
void AfterStep(double dt, StateVars& Q );

//void ConstructDiffTransformL( const StateVars& Q, const dTensorBC2& F,  dTensorBC2& Lstar, dTensorBC1& smax);

double GetCFL(double dt, double dtmax, const dTensorBC2& aux, const dTensorBC1& smax);

void BeforeFullTimeStep(double dt, const StateVars& Qold, StateVars& Qnew );
void AfterFullTimeStep (double dt, const StateVars& Qold, StateVars& Qnew );
// ------------------------------------------------------------

// ------------------------------------------------------------
// Solver that uses flux limiters
// ------------------------------------------------------------
void FinSolveUser( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC2& qnew = Qnew.ref_q  ();
    dTensorBC2&  aux = Qnew.ref_aux();

    // Time stepping information
    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number
    double t            = Qnew.get_t();     // Current time
    double dt           = dtv[1];           // Start with time step from last frame
    double cfl          = 0.0;              // current CFL number
    double dtmin        = dt;               // Counters for max and min time step taken
    double dtmax        = dt;

    // Grid information
    const int mx    = qnew.getsize(1);
    const int meqn  = qnew.getsize(2);
    const int maux  = aux.getsize(2);
    const int mbc   = qnew.getmbc();
    const int numel  = qnew.numel();

    // Maximum wave speed
    dTensorBC1    smax(mx, mbc);

    // Needed for rejecting a time step
    StateVars Qold( t, mx, meqn, maux, mbc );
    dTensorBC2& qold   = Qold.ref_q();
    dTensorBC2& auxold = Qold.ref_aux();

    // Allocate storage for this solver

    // Intermediate stages
    StateVars Qstar( t, mx, meqn, maux, mbc );
    dTensorBC2&   qstar = Qstar.ref_q();
    dTensorBC2& auxstar = Qstar.ref_aux();

    dTensorBC2  Lstar(mx, meqn, mbc); // right hand side (for Euler steps)
    dTensorBC2    F(mx, meqn, mbc );  // time-integrated flux

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
        if (n_step>nv)
        {
            cout << " Error in DogSolveUser.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold in case of rejecting a time step
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

            BeforeFullTimeStep(dt, Qold, Qnew);

            // ----------------------------------------------------------------
            // Take a single full time step of the solution
            // ----------------------------------------------------------------

            // Construct a "time-averaged" flux by considering Taylor
            // expansions at a single point.
            void ConstructDiffTransformL( double dt, StateVars& Q, dTensorBC1& smax, dTensorBC2& F, dTensorBC2& Lstar);
            ConstructDiffTransformL( dt, Qnew, smax, F, Lstar);

            // Now, check the fluxes to see if we need any limiting
            void LimitFluxes( const double dt, StateVars& Q, dTensorBC2& F);
            LimitFluxes( dt, Qnew, F);

            // Now use the time-averaged flux
            SetBndValues( Qnew );
            ConstructLxWL( aux, qnew, F, Lstar, smax); 

            // Update the solution:
#pragma omp parallel for
            for( int k=0; k < numel; k++ )
            {
                double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                qnew.vset(k, tmp );
            }
            Qnew.set_t( t );

            // ----------------------------------------------------------------

            // do any extra work      
            AfterFullTimeStep(dt, Qold, Qnew );

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

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
            if (cfl<=CFL_max) // accept
            { m_accept = 1; }
            else //reject
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

        }

        // compute (scalar) conservation values and print to file
        SetBndValues( Qnew );
        ConSoln     ( Qnew );

    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

}


void LimitFluxes( const double dt, StateVars& Q, dTensorBC2& F)
{

    SetBndValues( Q );
    dTensorBC2& q   = Q.ref_q();
    dTensorBC2& aux = Q.ref_aux();

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    // Multiresolution analysis parameter epsilon
    const double eps_mr  = 0.1;

    dTensorBC2 flag( mx, meqn, mbc );
    flag.setall(0.);

    // --------------------------------------------------------------------
    // Part I: Fill in values for the flag
    // --------------------------------------------------------------------
    int num_flagged_elems = 0;
    for( int i = 1-mbc+1; i <= mx+mbc-1; i++ )
    {

        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i+1,m) + q.get(i-1,m) );
            Qavg.set(m, tmp );

            if( abs( tmp - q.get(i,m) ) > eps_mr * dx )
            { 
                flag.set(i, m, 1); 
                num_flagged_elems += 1;     
                cout << "flagged elem i  = " << i << endl;
                F.set(i,m, q.get(i,m) );
            }

        }

    }
    // cout << "num_flagged_elems = " << num_flagged_elems << endl;

    // --------------------------------------------------------------------
    // Part II: Create buffer zone around each flagged element
    // --------------------------------------------------------------------
    dTensorBC2 buffed_flag( mx, meqn, mbc );
    buffed_flag.setall(0.);
    for( int i = 1; i <= mx; i++ )
    {

        for( int m=1; m <= meqn; m++ )
        {
            if( flag.get(i,m) > 0 )
            { 
                for( int j=-mbc; j <= mbc; j++ )
                {
                    buffed_flag.set(i+j, m, 1); 
                }
            }

        }

    }

    for( int i = 1-mbc; i <= mx+mbc; i++ )
    {

        for( int m=1; m <= meqn; m++ )
        {
            if( buffed_flag.get(i,m) > 0 && global_ini_params.get_mr_limiter() )
//          if( 1 )  // uncomment to force the "limiter" everywhere
            { 
//              cout << "limiting element i = " << i << endl;
                F.set(i,m, q.get(i,m) );
                F.set(i,m, q.get(i,m) - 0.5*dt );
            }

        }

    }

}
