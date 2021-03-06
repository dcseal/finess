#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "IniParams.h"

using namespace std;

// ------------------------------------------------------------
// Function definitions
void ConSoln( 
    const dTensorBC2& aux,
    const dTensorBC2& q, 
    double t, string outputdir);

// RK functions
void BeforeStep(double dt, dTensorBC2& aux, dTensorBC2& q);

void ConstructL(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        const dTensorBC2& F,  // <-- new term: integrated flux
        dTensorBC2& Lstar,
        dTensorBC1& smax);

void AfterStep(double dt, dTensorBC2& aux, dTensorBC2& q);

double GetCFL(double dt, double dtmax,
        const dTensorBC2& aux,
        const dTensorBC1& smax);

void BeforeFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);
void AfterFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);
// ------------------------------------------------------------

// ------------------------------------------------------------
// Taylor series integration
void ConstructIntegratedF( double dt, 
    dTensorBC2& aux, dTensorBC2& q,
    dTensorBC1& smax, dTensorBC2& F);

// ------------------------------------------------------------
// Multiderivative integration
void ConstructIntegratedF( double dt, 
    double alpha1, double beta1,
    dTensorBC2& aux1, dTensorBC2& q1,
    double alpha2, double beta2,
    dTensorBC2& aux2, dTensorBC2& q2,
    dTensorBC1& smax, dTensorBC2& F);

void ConstructIntegratedF( double dt, 
    double alpha1, double beta1, double charlie1,
    dTensorBC2& aux1, dTensorBC2& q1,
    double alpha2, double beta2, double charlie2,
    dTensorBC2& aux2, dTensorBC2& q2,
    dTensorBC1& smax, dTensorBC2& F);


void SetBndValues(dTensorBC2& aux, dTensorBC2& q);

// ------------------------------------------------------------
// Example function that describes how the main time stepping loop works.
// 
// This routine will terminate the code upon entrance.
// ------------------------------------------------------------
void DogSolveUser(
        dTensorBC2& aux, 
        dTensorBC2& qold,
        dTensorBC2& qnew,
        dTensorBC1& smax,
        double tstart, double tend,int nv, 
        double dtv[], const double cflv[],string outputdir)
{

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double CFL_max      = cflv[1];  // max   CFL number
    double CFL_target   = cflv[2];  // target CFL number
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const int mx    = qold.getsize(1);
    const int meqn  = qold.getsize(2);
    const int maux  = aux.getsize(2);
    const int mbc   = qnew.getmbc();

    // Total number of entries in the vector:
    const int numel = qold.numel();

    // Allocate storage for this solver
    dTensorBC2   qstar(mx, meqn, mbc);

    // extra aux arrays
    dTensorBC2  auxstar( mx, maux, mbc );

    dTensorBC2   Lstar(mx, meqn, mbc);

    // Time-averaged flux function
    dTensorBC2   F(mx, meqn, mbc);

    // Set initialize qstar and auxstar values
    qstar.copyfrom( qold );   
    auxstar.copyfrom( aux );

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step = 0;
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

        // copy qnew into qold
        qold.copyfrom( qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // ----------------------------------------------------------------
            BeforeFullTimeStep(dt, aux, aux, qold, qnew);


            SetBndValues(aux,      qnew);
            SetBndValues(auxstar, qstar);

            switch( global_ini_params.get_time_order() )
            {


                case 4:

                // -- Stage 1 -- //
    //          ConstructIntegratedF( 0.5*dt, aux, qnew, smax, F);

                ConstructIntegratedF( 0.5*dt, 
                    1.0, 0.5, aux,     qnew, 
                    0.0, 0.0, auxstar, qstar,
                    smax, F);
                ConstructL( aux, qnew, F, Lstar, smax);

                // Update the solution:
    #pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + 0.5*dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, auxstar, qstar );

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    1.0, (1.0/6.0), aux, qnew, 
                    0.0, (1.0/3.0), auxstar, qstar,
                    smax, F);
                ConstructL( auxstar, qstar, F, Lstar, smax);

                // Update the solution:
    #pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // Perform any extra work required:
                AfterStep(dt, auxstar, qstar );

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
                ConstructL( aux, qnew, F, Lstar, smax);

                // Update the solution:
    #pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + 2.0/5.0*dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }

                // Perform any extra work required:
                AfterStep(dt, auxstar, qstar );

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // -- Stage 2 -- //
                ConstructIntegratedF( dt, 
                    1.0, 0.5, 1./16.,         aux, qnew, 
                    0.0, 0.0, (5.0/48.0), auxstar, qstar,
                    smax, F);
                ConstructL( auxstar, qstar, F, Lstar, smax);

                // Update the solution:
    #pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }

                SetBndValues(aux,      qnew);
                SetBndValues(auxstar, qstar);

                // Perform any extra work required:
                AfterStep(dt, auxstar, qstar );

                break;

                default:
                printf("Error.  Time order %d not implemented for multiderivative\n", global_ini_params.get_time_order() );
                exit(1);

            }

            // do any extra work      
            AfterFullTimeStep(dt,aux,aux,qold,qnew);

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
                qnew.copyfrom( qold );
            }

        }

        // compute conservation and print to file
        ConSoln(aux, qnew, t, outputdir);

    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

}
