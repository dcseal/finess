#include "tensors.h"
#include "StateVars.h"
#include "IniParams.h"
#include <cmath>


using namespace std;


// Function prototypes
void applyBoundaryCondition( 
        const dTensorBC2& aux, 
        double alpha,   double beta, 
        double& A_n,    double& B_n );


void fastConvolve(
        dTensorBC2& aux, 
        double nu   );


// Function that is called before a full time step
// To avoid changing the basic structure of FINESS, we 
// will implement the implicit Wave equation solver here.
// In the actual time stepping loop, no changes are made 
// since the flux is set to zero.
//
// The implicit solver is detailed in
// [1]  Method of Lines Transpose: A Fast Implicit Wave Propagator
//      Matthew F. Causley, Andrew J. Christlieb, Yaman Guclu, Eric Wolf
//      arXiv:1306.6902
//      2013
//
void BeforeFullTimeStep(double dt, const StateVars& Qold, StateVars& Qnew )
{

    // Make a copy of the state variables
    Qnew.copyfrom(Qold);

    // The state variables are here 
    dTensorBC2& qvals   = Qnew.ref_q();
    // ... and the auxiliary variables are here
    dTensorBC2& aux     = Qnew.ref_aux();
    // We will be working with the auxiliary variables
    // in order to avoid changes to the basic program 
    // structure


    // Set some basic parameters
    int    n    = aux.getsize(1) - 1;   // In [1], the domain is 
                                        // discretized using n+1 points
    double xa   = aux.get(1,4);         // Domain boundary
    double xb   = aux.get(n+1,4);
    double c    = global_ini_params.get_cs_light(); // Speed of light
    double h    = global_ini_params.get_dx();       // Grid size



    // Some method-specific parameters (see [1] for details)
    // TODO: Move these to IniParams since they are used in several
    // places
    double beta     = 2.0;
    double nu       = beta*h/(c*dt);
    double alpha    = beta/(c*dt);

    // Some temporary storage
    int     ia      = 0;    // Loop inex
    double  u_nxt   = 0.0;  // Used to store solution at t^{n+1}
    // For applying boundary conditions
    double  A_n     = aux.get(0, 3);
    double  B_n     = aux.get(n+2, 3);


    // Implicit wave solver implementation
    // Step 1: 
    // Apply Green's function of modified Helmholtz operator 
    // (this is implemented in [1] as a fast convolution)
    fastConvolve(aux, nu);


    // Step 2: 
    // Apply boundary conitions
    // Note: only outflow boundary conditions available currently
    applyBoundaryCondition( aux, alpha, beta, A_n, B_n );
    // Store boundary constants for use at next time step
    aux.set(0, 3, A_n);
    aux.set(n+2, 3, B_n);


    // Step 3:
    // Compute solution at time t^{n+1}
    for (ia=1; ia<=n+1; ia++)  {
        // Solution at t^{n+1}
        u_nxt = -(pow(beta,2.0) - 2.0)*aux.get(ia, 1) - 
            aux.get(ia,2) + 
            (pow(beta,2.0)/2.0)*(   aux.get(ia,3) + 
                A_n*exp( -alpha*(aux.get(ia,4)-xa) ) + 
                B_n*exp( -alpha*(xb-aux.get(ia,4)) )    );
        // Update solution at t^{n-1}
        aux.set(ia, 2, aux.get(ia, 1) );
        aux.set(ia, 1, u_nxt);
    }


    // Copy solution to the state variables
    for( ia=1; ia<=n+1; ia++ )
        qvals.set( ia, 1, aux.get(ia, 1) );


    // Print out error in computed solution
    // Note: This is possible since we are simulating a manufactured
    // solution
    // u(x,t) = exp( -25 * (x-1/4+ct)^2 ) + exp( -25 * (x+1/4-ct)^2 )
    // Print out the inf. norm error
    double  error       =   0.0;
    double  error_max   =   0.0;
    double  t           =   Qnew.get_t() + dt;  // current time
    for( ia=1; ia<=n+1; ia++ )  {
        error = abs( 
                qvals.get(ia, 1) - 
                (   exp( -25.0*pow(aux.get(ia,4)+0.25-c*t, 2.0) ) + 
                    exp( -25.0*pow(aux.get(ia,4)-0.25+c*t, 2.0) )   )
                );
        if( error>error_max )
            error_max = error;
    }
    // cout << "Time step is " << dt << endl;
    cout << "Error in solution at time " << t << " is " 
            << error_max << endl;

}



// Function to impose (outflow) boundary conditions
// TODO: Other boundary conditions
// See [1] for details. This implements second order accurate outflow
// boundary conditions
void applyBoundaryCondition( 
        const dTensorBC2& aux, 
        const double alpha,   const double beta, 
        double& A_n,    double& B_n )   
{


    // Set some basic parameters
    int     n   = aux.getsize(1) - 1;   // Note: In [1], the domain 
                                        // divided into n+1 points
    double xa   = aux.get(1,4);         // Domain boundary
    double xb   = aux.get(n+1,4);

    // Temporary variables
    // For definitions of these, see [1]
    double  mu  = exp( -alpha*(xb-xa) );
    double  w_a = 0.0; 
    double  w_b = 0.0;
    double  gamma_0 =   ( 1.0 - exp(-beta) )/( pow(beta,2.0) ) - 
                        ( 1.0 + exp(-beta) )/( 2.0*beta );
    double  gamma_1 =   -2.0*( 1.0 - exp(-beta) )/( pow(beta,2.0) )
                            + (2.0/beta)*exp(-beta) + 1.0;
    double  gamma_2 =   ( 1.0 - exp(-beta) )/( pow(beta,2.0) ) + 
                        ( 1.0 - 3.0*exp(-beta) )/(2.0*beta) - 
                        exp(-beta);
    
    double  Gamma_0 =   0.5*gamma_0*pow(beta,2.0);
    double  Gamma_1 =   gamma_1 - gamma_0*( pow(beta,2.0) - 2.0 );
    double  Gamma_2 =   gamma_2 - gamma_0;


    // The boundary conditions are applied by solving a 2x2 system
    // These are the constants which define the right-hand side of the
    // system
    w_a     =   exp(-beta)*A_n + 
                Gamma_0*aux.get(1,3) +  // Gamma_0 * I[0]
                Gamma_1*aux.get(1,1) +  // Gamma_1 * u[0]
                Gamma_2*aux.get(1,2);   // Gamma_2 * u_prev[0]

    w_b     =   exp(-beta)*B_n + 
                Gamma_0*aux.get(n+1,3) +  // Gamma_0 * I[n]     
                Gamma_1*aux.get(n+1,1) +  // Gamma_1 * u[n]     
                Gamma_2*aux.get(n+1,2);   // Gamma_2 * u_prev[n]
                                        
    // This is where we solve the system
    A_n     =       ( (1.0-Gamma_0)*w_a + mu*Gamma_0*w_b ) / 
                    ( pow(1.0-Gamma_0,2.0) - pow(mu*Gamma_0,2.0) );
    B_n     =       ( (1.0-Gamma_0)*w_b + mu*Gamma_0*w_a ) / 
                    ( pow(1.0-Gamma_0,2.0) - pow(mu*Gamma_0,2.0) );                                

}



// The Maxwell solver requires the evaluation of a convolution
// operator (application of the Green's function of the modified
// Helmholtz operator)
// This can be evaluated efficiently in O(N) time (see [1] for details)
// Define the fast convolution operation here
void fastConvolve(
        dTensorBC2& aux, 
        const double nu   )  
{

    // Basic parameters
    int     ia  = 0;                    // Loop index
    int     n   = aux.getsize(1) - 1;   // In [1], the domain is 
                                        // discretized using n+1 points
    
    // Define some constants
    // These are used for evaluating the convolution integral
    // See [1] for details
    double  d   =   exp(-nu);
    double  P   =   1.0 - ( (1.0 - d)/nu );
    double  Q   =   -d + ( (1.0 - d)/nu );
    double  R   =   ( (1.0 - d)/pow(nu,2.0) ) - 
                    ( (1.0 + d)/(2.0*nu) ); 
    
    // Some temporary storage
    // TODO: I think this can be avoided; if not, move these to
    // auxiliary storage to avoid allocation at every function call
    double I_L[n+1], I_R[n+1], J_L[n+1], J_R[n+1];
    for( ia=0; ia<=n; ia++ )    {
        I_L[ia] = 0.0; I_R[ia] = 0.0; 
        J_L[ia] = 0.0; J_R[ia] = 0.0;
    }

    // The quadrature evaluation proceeds in three passes:
    // Left pass (see [1] for algorithm)
    I_L[0] = 0.0;
    for(ia=1; ia<=n; ia++)  {
        J_L[ia] = P*aux.get(ia+1,1) + Q*aux.get(ia,1);
        I_L[ia] = d*I_L[ia-1] + J_L[ia];
    }

    // Right pass
    I_R[0] = 0.0;
    for(ia=0; ia<=n-1; ia++)  {
        J_R[ia] = P*aux.get(ia+1,1) + Q*aux.get(ia+2,1);
        I_R[ia] = d*I_R[ia+1] + J_R[ia];
    }

    // Final pass
    aux.set(1,3,I_R[0]);
    for(ia=1; ia<=n-1; ia++)
        aux.set(ia+1,3,
                d*( I_L[ia-1] + I_R[ia+1] ) + 
                2.0*P*aux.get(ia+1,1) + 
                Q*( aux.get(ia+2,1) + aux.get(ia,1) ) + 
                2.0*R*( aux.get(ia+2,1) - 2.0*aux.get(ia+1,1) + 
                        aux.get(ia,1) )     
               );
    aux.set(n+1,3,I_L[n]);


}


