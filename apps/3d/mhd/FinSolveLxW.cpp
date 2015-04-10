#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "FinSolveLxW.h"     // functions directly called from this function
#include "IniParams.h"
#include "ConstructHJ_L.h"

using namespace std;

void FinSolveLxW ( StateVars& Qnew, double tend, double dtv[] )
{
    void (*ConstructHJ_L)(const StateVars& Q, dTensorBC4& Lauxstar, double dt) = ConstructHJ_L_Order3;

    dTensorBC4& qnew = Qnew.ref_q  ();
    dTensorBC4&  aux = Qnew.ref_aux();

    // Time stepping information
    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number
    double t                  = Qnew.get_t();
    double dt                 = dtv[1];   // Start with time step from last frame
    double cfl                = 0.0;      // current CFL number
    double dtmin              = dt;       // Counters for max and min time step taken
    double dtmax              = dt;

    // Grid information
    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();
    const int mz   = global_ini_params.get_mz();
    const int meqn = global_ini_params.get_meqn();
    const int maux = global_ini_params.get_maux();
    const int mbc  = global_ini_params.get_mbc();
    const int numel = qnew.numel();
    const int numel_aux = aux.numel();

    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();
    const double dz = global_ini_params.get_dz();

    // Maximum wave speed at each flux interface value.  Note that the size of
    // this tensor is one longer in each direction.
    dTensorBC4 smax( mx+1, my+1, mz+1, 3, mbc );           

    // Needed for rejecting a time step
    StateVars Qold( t, mx, my, mz, meqn, maux, mbc );
    dTensorBC4& qold   = Qold.ref_q();
    dTensorBC4& auxold = Qold.ref_aux();
    Qold.copyfrom( Qnew );

    // Allocate storage for this solver
    dTensorBC4       F(mx, my, mz, meqn, mbc );  // time-integrated flux
    dTensorBC4       G(mx, my, mz, meqn, mbc );  // time-integrated flux
    dTensorBC4       H(mx, my, mz, meqn, mbc );  // time-integrated flux
    dTensorBC4   Lstar(mx, my, mz, meqn, mbc);   // Right hand side of ODE
    dTensorBC4   Lauxstar(mx, my, mz, maux, mbc);

    // Storage for the MPP limiter
    dTensorBC4* fhat;
    dTensorBC4* fLF;
    dTensorBC4* ghat;
    dTensorBC4* gLF;
    dTensorBC4* hhat;
    dTensorBC4* hLF;
    if( global_ini_params.get_mpp_limiter() )
    {
        fhat = new dTensorBC4( mx+1, my, mz, meqn, mbc );
        fLF  = new dTensorBC4( mx+1, my, mz, meqn, mbc );

        ghat = new dTensorBC4( mx, my+1, mz, meqn, mbc );
        gLF  = new dTensorBC4( mx, my+1, mz, meqn, mbc );

        hhat = new dTensorBC4( mx, my, mz+1, meqn, mbc );
        hLF  = new dTensorBC4( mx, my, mz+1, meqn, mbc );
    }


    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;                           // Number of time steps taken
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

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

            // ---------------------------------------------------------
            // Take a full time step of size dt
            BeforeStep(dt, Qnew);
            SetBndValues(Qnew);
            ConstructIntegratedR( dt, Qnew, smax, F, G, H);


            if( global_ini_params.get_mpp_limiter() )
            {

                // Construct the high-order flux
                ConstructLxWL( Qnew, F, G, H, fhat, ghat, hhat, Lstar, smax );

                 if(global_ini_params.get_constrained_transport() == true)
                    ConstructHJ_L(Qnew, Lauxstar, dt);
               // Construct the low-order flux
                ConstructLFL( dt, Qnew, *fLF, *gLF, *hLF, Lstar, smax );

                // Limit the high-order flux
                ApplyMPPLimiter3D( dt, qnew, *fLF, *gLF, *hLF, *fhat, *ghat, *hhat);

                // Update the solution:
#pragma omp parallel for
                for( int i=1; i <= mx; i++   )
                    for( int j=1; j <= my; j++   )
                        for( int k=1; k <= mz; k++   ){
                            for( int m=1; m <= meqn; m++ )
                            {
                                double tmp = fhat->get(i+1,j,k,m)-fhat->get(i,j,k,m) ;
                                qnew.set(i, j, k, m, qnew.get(i,j,k,m) - (dt/dx)*tmp );
                                tmp = ghat->get(i,j+1,k,m)-ghat->get(i,j,k,m);
                                qnew.set(i, j, k, m, qnew.get(i,j,k,m) - (dt/dy)*tmp );
                                tmp = hhat->get(i,j,k+1,m)-hhat->get(i,j,k,m) ;
                                qnew.set(i, j, k, m, qnew.get(i,j,k,m) - (dt/dz)*tmp );

                                // Test the LF solver by uncommenting this chunk
//                                double tmp = Lstar.get(i,j,k,m);
//                                qnew.set(i, j, k, m, qnew.get(i,j,k,m) + dt*tmp );

                            }
                            double rho = qnew.get(i, j, k, 1);
                            double u1 = qnew.get(i,j,k,2)/rho;
                            double u2 = qnew.get(i,j,k,3)/rho;
                            double u3 = qnew.get(i,j,k,4)/rho;
                            double energy = qnew.get(i,j,k,5);
                            double B1 = qnew.get(i,j,k,6);
                            double B2 = qnew.get(i,j,k,7);
                            double B3 = qnew.get(i,j,k,8);
                            double gamma = global_ini_params.get_gamma();
                            double p = (gamma - 1.0) * (energy - rho*(u1*u1+u2*u2+u3*u3)*0.5 - (B1*B1+B2*B2+B3*B3)*0.5);
                            //			if(p < 0.0)
                            //			    cerr << "( " << i << ", " << j << ", " << k << ")  p=" << p << "\n";
                        }
#pragma omp parallel for
                for(int k =  0; k < numel_aux; k++){
                    double tmp = aux.vget(k) + dt*Lauxstar.vget(k);
                    aux.vset(k, tmp);
                }

            }
            else
            { 

                // Construct RHS
                ConstructLxWL( Qnew, F, G, H, Lstar, smax);

                if(global_ini_params.get_constrained_transport() == true)
                    ConstructHJ_L(Qnew, Lauxstar, dt);


                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
#pragma omp parallel for
                for(int k =  0; k < numel_aux; k++){
                    double tmp = aux.vget(k) + dt*Lauxstar.vget(k);
                    aux.vset(k, tmp);
                }

            }

            Qnew.set_t( Qnew.get_t() + dt );

            // Perform any extra work required:
            AfterStep(dt, Qnew);
            // ---------------------------------------------------------

            // do any extra work
            AfterFullTimeStep(dt, Qold, Qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveLxW3D ... Step" << setw(5) << n_step;
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
                    cout<<"FinSolveLxW3D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                Qnew.copyfrom( Qold  );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln(Qnew);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    // Clean up allocated memory
    if( global_ini_params.get_mpp_limiter() )
    {
        delete fhat;
        delete fLF;
        delete ghat;
        delete gLF;
        delete hhat;
        delete hLF;
    }
}
