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


// TODO - move this into a separate file.  See note below
void ConstructDiffTransformL( const double dt, StateVars& Q, dTensorBC1& smax, dTensorBC2& F);

double GetCFL(double dt, double dtmax, const dTensorBC2& aux, const dTensorBC1& smax);

void BeforeFullTimeStep(double dt, const StateVars& Qold, StateVars& Qnew );
void AfterFullTimeStep (double dt, const StateVars& Qold, StateVars& Qnew );
// ------------------------------------------------------------

// ------------------------------------------------------------
// Example function that describes how the main time stepping loop works.
// 
// This routine will terminate the code upon entrance.
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
            ConstructDiffTransformL( dt, Qnew, smax, F);

            // Now use the time-averaged flux
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

// TODO - move this into a separate file, ConstructDiffTransformL.cpp and
// correctly link this.  When I try compiling this with the following lines
// in lib/Makefile.defs:
//
//  # place to add objects and sources
//  #
//  APP_LIB_OBJECTS = \
//    $(APP_LIB)/ConstructDiffTransformL.o
//
//  #
//  APP_LIB_SOURCES = \
//    $(APP_LIB)/ConstructDiffTransformL.cpp
//
//  #
//
// I end up with two references to ConstructDiffTransformL.o in the linking
// step.  This seems to work in other parts of the code though ...

#include "CentralDifferences.h"

// Time expanded state variable, q using discrete transform.
// See: Multi-moment ADER-Taylor methods for systems of conservation laws
// with source terms in one dimension 
//
// Euler equations.
//
// See also: $FINESS/lib/1d/ConstructIntegratedF.cpp.
void ConstructDiffTransformL( const double dt, StateVars& Q, dTensorBC1& smax, dTensorBC2& F)
{

    dTensorBC2& q   = Q.ref_q();
    dTensorBC2& aux = Q.ref_aux();

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    const double gamma = global_ini_params.get_gamma();


    const int mpts_sten       = global_ini_params.get_space_order(); assert_eq( mpts_sten%2, 1 );

    const int mbc_small       = (mbc+1)/2;
    const int half_mpts_sten  = (mpts_sten+1)/2;          // assert_eq( half_mpts_sten, 3 );
    const int MAX_DERIVS      = mpts_sten;
    const int MAX_FLUX_DERIVS = mpts_sten-1;

    // Quadrature rules for numerically evaluating the integral of the flux
    const int NumQuad = 9;
    dTensor1 w1d( NumQuad );
    dTensor1 x1d( NumQuad );

    void setGaussLobattoPoints1d( dTensor1& w1d, dTensor1& x1d);
    void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);


    //printf("Here=%e \n",gamma);


    //if( mpts_sten < 9 )
    if( NumQuad < 9 )
        setGaussLobattoPoints1d( w1d, x1d);
    else
        setGaussPoints1d( w1d, x1d );

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    {

        // Physical location for this current value:
        dTensor1 xpts( 1 );
        xpts.set( 1, xlow + double(i)*dx - 0.5*dx );

        // Save the flux function:
        dTensor2 qvals( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                qvals.set( m, r, q.get( i+s, m ) );
                s++;
            }
        }

        // Storage for local derivatives:
        dTensor2 qderivs ( meqn, MAX_DERIVS );

        // Compute a FD approximation to the derivatives:
        // Central difference routine (depends on spatial order!)
        void (*CentralDifferences)( double dx, const dTensor2& f, dTensor2& fderivs) = GetCentralDifferences();
        CentralDifferences( dx, qvals, qderivs );

        dTensor3 Gs( 1, MAX_DERIVS, MAX_DERIVS );
        dTensor3 R( 1, MAX_DERIVS, MAX_DERIVS );
        dTensor3 Gr( 2, MAX_DERIVS, MAX_DERIVS );
        dTensor3 P( 1, MAX_DERIVS, MAX_DERIVS );

        // Initialize all arrays to zero (maybe done already...i should check the constructor)
        Gs.setall( 0. );
        R.setall ( 0. );
        Gr.setall( 0. );
        P.setall ( 0. );


//      for( int h=1; h <= MAX_DERIVS; h++ )
//      for( int k=1; k <= MAX_DERIVS; k++ )
//      {
//          for( int me =1; me <= meqn; me++ )
//          {
//              Q_mixed_derivs.set( me, h, k, 0.0 );
//          }
//          Gs.set(1,h,k,0.0);
//          R.set(1,h,k,0.0);
//          Gr.set(1,h,k,0.0);
//          Gr.set(2,h,k,0.0);
//          Ge.set(1,h,k,0.0);
//          P.set(1,h,k,0.0);
//      }

        // Save all of the "zeroth" time derivatives.
        dTensor3 Q_mixed_derivs( meqn, MAX_DERIVS, MAX_DERIVS );
        for( int me =1; me <= meqn; me++ )
        {
            for( int h=1; h <= MAX_DERIVS; h++ )
            { 
                Q_mixed_derivs.set( me, h, 1, qderivs.get(me, h)/factorial[h-1] ); 
            }
        }

        //Need to deal with k=0 specially because of the 1/rho term...
        for(int k=0;k<1;k++)
        {

          //Populate time derivatives of (rho u)^2...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 2, r+1, s+1)*Q_mixed_derivs.get(2, h+1-r, k+1-s );
                }
                Gs.set( 1, h+1, k+1, tmp );
          }

          R.set( 1, 1, 1, 1.0/(Q_mixed_derivs.get(1,1,1)) );

          //Populate time derivatives of 1/rho...
          for( int h=1; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                R.set(1,h+1,k+1,0.0);
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1));
                R.set( 1, h+1, k+1, tmp );
          }

          //populate the time derivatives of (rho u)^2/rho and u
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Gs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                      tmp1 += Q_mixed_derivs.get( 2, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                Gr.set( 1, h+1, k+1, tmp );
                Gr.set( 2, h+1, k+1, tmp1 );
          }  

          //Populate the time derivatives of pressure
          
         
          
          //Populate time derivatives of pressure
          /*for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                tmp = 0.4*(Q_mixed_derivs.get( 5, h+1, k+1)-0.5*Gr.get(1,h+1,k+1));
                P.set( 1, h+1, k+1, tmp );
          }*/
          
          //Populate time derivatives of pressure
         
          P.set(1,1,1, pow(Q_mixed_derivs.get(1,1,1),gamma));
           

          for( int h=0; h < MAX_DERIVS-1; h++ )
          {
                //double tmp = Q_mixed_derivs.get(1,h+1,k+1);
                //for( int r=0; r <= h; r++ )
                //for( int s=0; s <= k; s++ )
                //{
                //      tmp += Q_mixed_;
                //}
                double tmp = 0.0;
                double tmp1 = 0.0;
              
                for( int r=1; r <= h; r++ )
                for( int s=1; s <= k; s++ )
                {
                      tmp += (h+1-r)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+2,k-s+1);
                      tmp1 += (h+1-r)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+2,k-s+1);
                }

                for( int r=1; r <= h; r++ )
                for( int s=0; s <= 0; s++ )
                {     
                      tmp += (h+1-r)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+2,k-s+1);
                      tmp1 += (h+1-r)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+2,k-s+1);
                }

                for( int r=0; r <= 0; r++ )
                for( int s=1; s <= k; s++ )
                {     
                      tmp += (h+1-r)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+2,k-s+1);
                      tmp1 += (h+1-r)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+2,k-s+1);
                }

                tmp+=(h+1)*P.get(1,1,1)*Q_mixed_derivs.get(1,h+2,k+1);
                P.set( 1, h+2, k+1, 1.0/((h+1)*Q_mixed_derivs.get(1,1,1))*(gamma*tmp-tmp1) );
          }

          
          for( int h=0; h < MAX_FLUX_DERIVS-k; h++ )
          {

           //compute conserved quantity derivatives...
           Q_mixed_derivs.set(1,h+1,k+2,-((double)h+1.0)/(k+1)*Q_mixed_derivs.get(2,h+2,k+1));
           Q_mixed_derivs.set(2,h+1,k+2,-((double)h+1.0)/(k+1)*(Gr.get(1,h+2,k+1)+P.get(1,h+2,k+1)));
           }

       }


        // Recursive relationship goes here! Compute higher order derivatives //
        for( int k=1; k < MAX_FLUX_DERIVS;   k++ )      
        {
        //Populate time derivatives of (rho u)^2...
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 2, r+1, s+1)*Q_mixed_derivs.get(2, h+1-r, k+1-s );
                }
                Gs.set( 1, h+1, k+1, tmp );
        }
        //Populate time derivatives of 1/rho...
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp = 0.0;
                R.set(1,h+1,k+1,0.0);
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1));
                R.set( 1, h+1, k+1, tmp );
        }

        //populate time derivatives of (rho u)^2/rho and u
        for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp = 0.0;
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Gs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                      tmp1 += Q_mixed_derivs.get( 2, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                Gr.set( 1, h+1, k+1, tmp );
                Gr.set( 2, h+1, k+1, tmp1 );
        }
        //populate time derivatives of P
        /*for( int h=0; h < MAX_DERIVS; h++ )
        {
                double tmp = 0.0;
                tmp = Q_mixed_derivs.get(1,h+1,k+1);
                P.set( 1, h+1, k+1, tmp );
        }*/


          int km1=k-1;
          for( int h=0; h < MAX_DERIVS-1; h++ )
          {     
                //double tmp = Q_mixed_derivs.get(1,h+1,k+1);
                //for( int r=0; r <= h; r++ )
                //for( int s=0; s <= k; s++ )
                //{
                //      tmp += Q_mixed_;
                //}
                double tmp = 0.0;
                double tmp1 = 0.0;
                
                for( int r=1; r <= h; r++ )
                for( int s=1; s <= km1; s++ )
                {     
                      tmp += (km1+1-s)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+1,km1-s+2);
                      tmp1 += (km1+1-s)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+1,km1-s+2);
                }
                
                for( int r=1; r <= h; r++ )
                for( int s=0; s <= 0; s++ )
                {     
                      tmp += (km1+1-s)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+1,km1-s+2);
                      tmp1 += (km1+1-s)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+1,km1-s+2);
                }
                
                for( int r=0; r <= 0; r++ )
                for( int s=1; s <= km1; s++ )
                {     
                      tmp += (km1+1-s)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+1,km1-s+2);
                      tmp1 += (km1+1-s)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+1,km1-s+2);
                }
                
                tmp+=(km1+1)*P.get(1,1,1)*Q_mixed_derivs.get(1,h+1,km1+2);
                P.set( 1, h+1, km1+2, 1.0/((km1+1)*Q_mixed_derivs.get(1,1,1))*(gamma*tmp-tmp1) );
          }


         /*
          for( int h=0; h < MAX_DERIVS-1; h++ )
          {
                //double tmp = Q_mixed_derivs.get(1,h+1,k+1);
                //for( int r=0; r <= h; r++ )
                //for( int s=0; s <= k; s++ )
                //{
                //      tmp += Q_mixed_;
                //}
                double tmp = 0.0;
                double tmp1 = 0.0;

                for( int r=1; r <= h; r++ )
                for( int s=1; s <= k; s++ )
                {     
                      tmp += (h+1-r)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+2,k-s+1);
                      tmp1 += (h+1-r)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+2,k-s+1);
                }

                for( int r=1; r <= h; r++ )
                for( int s=0; s <= 0; s++ )
                {
                      tmp += (h+1-r)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+2,k-s+1);
                      tmp1 += (h+1-r)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+2,k-s+1);
                }

                for( int r=0; r <= 0; r++ )
                for( int s=1; s <= k; s++ )
                {
                      tmp += (h+1-r)*P.get(1,r+1,s+1)*Q_mixed_derivs.get(1,h-r+2,k-s+1);
                      tmp1 += (h+1-r)*Q_mixed_derivs.get(1,r+1,s+1)*P.get(1,h-r+2,k-s+1);
                }

                tmp+=(h+1)*P.get(1,1,1)*Q_mixed_derivs.get(1,h+2,k+1);
                P.set( 1, h+2, k+1, 1.0/((h+1)*Q_mixed_derivs.get(1,1,1))*(gamma*tmp-tmp1) );
          }*/


          for( int h=0; h < MAX_FLUX_DERIVS-k; h++ )
          {

           //compute conserved quantity derivatives...
           Q_mixed_derivs.set(1,h+1,k+2,-((double)h+1.0)/(k+1)*Q_mixed_derivs.get(2,h+2,k+1));
           Q_mixed_derivs.set(2,h+1,k+2,-((double)h+1.0)/(k+1)*(Gr.get(1,h+2,k+1)+P.get(1,h+2,k+1)));
           }
         }

        //Compute the highest time-order terms of the auxillary quantities because we will need these for the Euler fluxes 
        for(int k=MAX_FLUX_DERIVS-1;k<MAX_DERIVS;k++)
        {

          //Populate time derivatives of (rho u)^2...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 2, r+1, s+1)*Q_mixed_derivs.get(2, h+1-r, k+1-s );
                }
                Gs.set( 1, h+1, k+1, tmp );
          }

          //Populate time derivatives of 1/rho...
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                R.set(1,h+1,k+1,0.0);
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Q_mixed_derivs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                tmp *=-1.0/(Q_mixed_derivs.get(1,1,1));
                R.set( 1, h+1, k+1, tmp );
          }

          //populate time derivatives of (rho u)^2/rho and u
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                double tmp1 = 0.0;
                for( int r=0; r <= h; r++ )
                for( int s=0; s <= k; s++ )
                {
                      tmp += Gs.get( 1, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                      tmp1 += Q_mixed_derivs.get( 2, r+1, s+1)*R.get(1, h+1-r, k+1-s );
                }
                Gr.set( 1, h+1, k+1, tmp );
                Gr.set( 2, h+1, k+1, tmp1 );
          }  

          //populate time derivatives of pressure
          for( int h=0; h < MAX_DERIVS; h++ )
          {
                double tmp = 0.0;
                tmp = Q_mixed_derivs.get(1,h+1,k+1);
                P.set( 1, h+1, k+1, tmp );
          }
       

       }

        // Construct the time-averaged flux.
        double rho_ta_flux  = 0.;
        double u1_ta_flux   = 0.;
        double E_ta_flux    = 0.;
        for( int mq=1; mq <= NumQuad; mq++ )
        {
            // Evaluate q at this quadrature point
            double rho  = Q_mixed_derivs.get(1,1,1);
            double ru1  = Q_mixed_derivs.get(2,1,1);
            double P1   = P.get(1,1,1);
            for( int k=1; k < MAX_DERIVS; k++ )
            { 
                rho += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 1, 1, k+1 ); 
                ru1 += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*Q_mixed_derivs.get( 2, 1, k+1 ); 

                //the pressure function 
                P1 += ( pow(0.5*dt*( 1.0 + x1d.get(mq) ), k) )*P.get( 1, 1, k+1 ); 
            }
            rho_ta_flux += (0.5*w1d.get( mq )) * ( ru1 );
            u1_ta_flux  += (0.5*w1d.get( mq )) * ( ru1*ru1/rho+P1);
        }

        F.set( i, 1, rho_ta_flux ); // flux for density
        F.set( i, 2, u1_ta_flux );  // flux for 1-component of momentum

    }

}
