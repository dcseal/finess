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
// this routine supports a third and fourth-order two derivative, and a 
// fifth-order three-derivative method.
//
// One needs to have a Jacobian defined in DFlux, as well as the Hessian in
// D2FluxFunc (for the three-derivative method) in order to be able to use 
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

    // Intermediate storage (for three stage methods)
    StateVars Q2star( t, mx, meqn, maux, mbc );
    dTensorBC2& q2star   = Q2star.ref_q();
    dTensorBC2& aux2star = Q2star.ref_aux();
    Q2star.copyfrom( Qnew );

    // Right hand side of ODE and time-averaged flux function
    dTensorBC2   Lstar(mx, meqn, mbc);
    dTensorBC2   F(mx, meqn, mbc);

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

                case 3:

                // In this section, we have two third-order methods with and
                // without Shu-Osher decomposition:
                //
                // (1) (Without a “Shu-Osher” decomposition) 
                //
                // y* = y^n + 2/3*dt W( f^n + dt/3 \dot{f}^n )
                // u^{n+1} = u^n + dt * W( 5/8 f^n + 3/8 f* + dt/8 *( \dof{f} + \dot{f*} ) )
                //
                // This First method exhibits an increase in total variation
                // for any positive CFL number.
                //
                // (2) (With a “Shu-Osher” decomposition)
                //
                // y* = y^n + 2/3*dt W( f^n + dt/3 \dot{f}^n )
                // u^{n+1} = 7/16 u^n + 9/16 y* + dt * W( 1/4 f^n + 3/8 f* + dt/8 \dof{f*} )


                // ----------------------------- //
                // No Shu-Osher decomposition
                // ----------------------------- //

                // -- Stage 1 -- //
/*
                ConstructIntegratedF( dt, 
                    A21, Ahat21, Qnew, 
                    0.0, 0.0,    Qstar,
                    smax, F);

                // Update the solution:
                ConstructLxWL( aux, qnew, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + (2./3.)*dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //

                // Construct a new right hand side
                ConstructIntegratedF( dt, 
                    b1, bhat1, Qnew, 
                    b2, bhat2, Qstar,
                    smax, F);
                ConstructLxWL( auxstar, qstar, F, Lstar, smax);

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
*/

                // ----------------------------- //
                // With a Shu-Osher decomposition
                // ----------------------------- //

                // -- Stage 1 -- //

                ConstructIntegratedF( dt, 
                    1.0, 0.5, Qnew, 
                    0.0, 0.0, Qstar,
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

                // Update the solution:
                ConstructIntegratedF( dt, 
                    0.25,      0.0, aux, qnew, 
                    3.0/8.0, 1./8., auxstar, qstar,
                    smax, F);

                ConstructLxWL( auxstar, qstar, F, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = (7.*qold.vget( k ) + 9.*qstar.vget(k))/16.;
                    tmp += dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qnew.get_t() + dt );

                SetBndValues( Qnew  );
                SetBndValues( Qstar );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                break;

/*

                // Update the solution:
                ConstructL( Qnew, Lstar, smax);
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

                // This is a three-stage, fifth-order non-SSP method.  It has 
                // negative coefficients. Coefficients can be found
                // in D.C. Seal, Y. Güçlü and A.J. Christlieb. "High-order
                // multiderivative time integrators for hyperbolic
                // conservation laws", J. Sci. Comp., Vol. 60, Issue 1, pp
                // 101-140, 2014.

                // -- Stage 1 -- //
                ConstructIntegratedF( dt, 
                    1.0, 0.25, Qnew, 
                    0.0, 0.0,  Qstar,
                    smax, F);
                ConstructLxWL( aux, qnew, F, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + (2.0/5.0)*dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }
                Qstar.set_t( Qnew.get_t() + (2.0/5.0)*dt );

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues( Qnew   );
                SetBndValues( Qstar  );
                SetBndValues( Q2star );

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    1.0, (1.0/6.0), Qnew, 
                    0.0, (1.0/3.0), Qstar,
                    smax, F);
                ConstructLxWL( aux, qnew, F, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    q2star.vset(k, tmp );
                }
                Q2star.set_t( Qold.get_t() + dt );

                SetBndValues( Qnew   );
                SetBndValues( Qstar  );
                SetBndValues( Q2star );

                // Perform any extra work required:
                AfterStep(dt, Q2star );

                // -- Stage 3 (final step) -- //
                ConstructIntegratedF( dt, 
                    1.0,  1.0/8.0,     Qnew, 
                    0.0,  25./72.0,   Qstar,
                    0.0,  1.0/36.0,  Q2star,
                    smax, F);
                ConstructLxWL( aux, qnew, F, Lstar, smax);


                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qold.get_t() + dt );

                SetBndValues( Qnew   );
                SetBndValues( Qstar  );
                SetBndValues( Q2star );

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                /*
                 * TODO - this section should be removed: it was the wrong
                 * thing to try and optimize.  -DS
                 */

                /*
                // Coeffients chosen to optimize region of absolute stability along the
                // imaginary axis.
                //
                // rho = 8.209945182837015e-02 chosen to maximize range of abs. stab. region

                // -- Stage 1 -- //
                ConstructIntegratedF( 2.0/5.0*dt, 
                    1.0, 0.5, 125./8.*8.209945182837015e-02, Qnew, 
                    0.0, 0.0, 0.0,                           Qstar,
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
                    1.0, 0.5, (1.0/16.0), Qnew, 
                    0.0, 0.0, (5.0/48.0), Qstar,
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
