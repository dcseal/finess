#include "tensors.h"
#include "StateVars.h"
#include <iostream>
#include <cmath>
#include "IniParams.h"
#include "Components.h"


using namespace std;




// The Maxwell solver requires the evaluation of a convolution
// operator (application of the Green's function of the modified
// Helmholtz operator)
// This can be evaluated efficiently in O(N) time (see [1] for details)
// Define the fast convolution operation here
// Implementation is identical to Matt's Matlab function, but restricted 
// to 2nd order accuracy
// TODO: Generalize to arbitrary order

// Here is the documentation from Matt's Matlab code
/*
% 
% ______________________________Description______________________________ %
% Fast convolution for the integral
% 
%           I[u](x) = \alpha \int_a^b  u(y) e^{-\alpha|x-y|}dy,
% 
% using N+1 arbitrarily spaced points, and quadrature of order M.
% 
% ________________________________Inputs_________________________________ %
% u   a vector of length N+1, evaluated at N+1 arbitrarily spaced x points
% nu  a vector of length N, where nu_j = alpha*(x_{j} - x_{j-1})
% 
% ________________________________Outputs________________________________ %
% I   a vector of length N+1, evaluated at N+1 arbitrarily spaced x points
% wL  an array of size (N+1)x(M+1) containing weights for IL
% wR  an array of size (N+1)x(M+1) containing weights for IR
% 
% __________________________Problem Description__________________________ %
% Let the collocation points be labelled as {a=x_0,x_1,...,x_N,x_N=b},
% where h_j =x_j-x_{j-1} may vary depending on j. Now, we decompose the
% integral into I = IL + IR, such that
% 
%  IL_j     = d_j IL_{j-1} + nu_j int_0^1 u(x_{j  }-z h_j)e^{-nu_j z}dz,
%  IR_{j-1} = d_j IR_{j  } + nu_j int_0^1 u(x_{j-1}+z h_j)e^{-nu_j z}dz,
% 
% where IL_0 = IR_N = 0, and d_j = e^{-nu_j}. These equations are exact.
% 
% _________________________Algorithm Description_________________________ %
% Now, upon replacing u with a Lagrange interpolant of degree M,
% 
% p(z) = sum_{m=-M/2}^{M/2}  p_m(z) u_{j+m}.
% 
% Near the endpoints, the stencil is adjusted accordingly. Now, we
% integrate the expression analytically, and get
% 
% IL_j = d_j IL_{j-1} +sum_{m=-M/2}^{M/2} wL_{j,m} u_{j+m}
% IR_j = d_j IR_{j+1} +sum_{m=-M/2}^{M/2} wR_{j,m} u_{j+m},
% 
% _______________________________________________________________________ %
*/

void fastConvolve(
        double I[], 
        const double u[],
        const int N,  
        const double nu)  
{


    // Basic parameters
    const int   M   = 2;                    // Order of approximation
    double      h   = global_ini_params.get_dx();       // Grid size 
    
    // Temporary variables
    int     k   = 0;                    // Loop index
    int     m   = 0;                    // Loop index
    double  IL  =   0.0;
    double  IR  =   0.0;
    int     M2  =   (int)(M/2);
    
    double  iL[M+1];                    // Initialize iL
    for( k=0; k<=M; k++ )
        iL[k] = 1.0 + (double) k;

    double  iR[M+1];                    // Initialize iR
    for( k=-M; k<=0; k++ )
        iR[k+M] = (double) N + 1.0 + (double) k;

    double  iC[2*M2+1];                 // Initialize iR
    for( k=-M2; k<=M2; k++ )
        iC[k+M2] = k;

    // Initialize integral solution, I
    for( k=0; k<=N; k++ )
        I[k] = 0.0;

    // Initialize Weighted nodes
    double nu_vec[N], d[N];
    for( k=0; k<N; k++ )    {
        nu_vec[k]   = nu;
        d[k]        = exp(-nu);
    }


    // Code to generate exponential weights
    // Here is Matt's documentation from his Matlab code
    /*
    % 
    % ______________________________Description______________________________ %
    % Construct the local weights wL and wR for integrating the local
    % contributions to the integrals
    % 
    %  JL_j     = nu_j int_0^1 u(x_{j  }-z h_j)e^{-nu_j z}dz,
    %  JR_{j-1} = nu_j int_0^1 u(x_{j-1}+z h_j)e^{-nu_j z}dz,
    % 
    % up to order M.
    % ________________________________Inputs_________________________________ %
    % nu  a vector of length N, where nu_j = alpha*(x_{j} - x_{j-1})
    % M   an even integer, which is the order of the integration
    % 
    % ________________________________Outputs________________________________ %
    % wL  an array of size [N,M+1] storing the weights of JL for each j
    % wR  an array of size [N,M+1] storing the weights of JR for each j
    % 
    % __________________________Problem Description__________________________ %
    % Upon discretization, u(x_j-z h_j) is replaced by the polynomial of degree
    % M. The integrals are then approximated with quadrature, and the weights
    % can be precomputed, so that
    % 
    %  JL_j     = sum_{m=-M/2}^{M/2} wL_{jm}u_{j+m} 
    %  JR_j     = sum_{m=-M/2}^{M/2} wR_{jm}u_{j+m}
    % 
    % _________________________Algorithm Description_________________________ %
    % Using nu_{j+m}, the Vandermonde matrix is constructed locally, and it's
    % inverse is multiplied by the coefficients phi_m, which are defined by
    % 
    % 
    % phi_m(nu) = nu int_0^1 z^m exp(-nu z)dz
    % 
    % This matrix vector product precisely defines the weights wL and wR.
    %
    % _______________________________________________________________________ %
    */
    double  wL[N][M+1], wR[N][M+1];
    // First, get exponential coefficients
    double  phi[M+1];
    if (nu<40.0) {
        double t = 1.0;
        double s = 0.0;
               k = 1;

        while (t>1.0e-15)   {
            t = t*nu/( (double) M + (double) k );
            s = s+t;
            k = k+1;
        }
        phi[M] = exp(-nu)*s;
    }
    else    {
        double t = 1.0;
        double s = 0.0;

        for(k=1; k<=M; k++) {
            t = t*nu / (double) k;
            s = s+t;
        }
        phi[M] = 2.0 / (nu*nu)*(1.0-exp(-nu)*s);
    }

    for( m=M; m>=1; m-- )
        phi[m-1] = nu / (double)m * ( phi[m] + exp(-nu) );


    // Now, construct weights
    // wL
    for( k=0; k<N-1; k++ )    {
        wL[k][0] = 0.0 * phi[0] + 0.5 * phi[1] + 0.5 * phi[2];
        wL[k][1] = 1.0 * phi[0] + 0.0 * phi[1] - 1.0 * phi[2];
        wL[k][2] = 0.0 * phi[0] - 0.5 * phi[1] + 0.5 * phi[2];
    }

    wL[N-1][0] = 0.0 * phi[0] - 0.5 * phi[1] + 0.5 * phi[2];
    wL[N-1][1] = 0.0 * phi[0] + 2.0 * phi[1] - 1.0 * phi[2];
    wL[N-1][2] = 1.0 * phi[0] - 1.5 * phi[1] + 0.5 * phi[2];

    // wR
    for( k=1; k<N; k++ )    {
        wR[k][0] = 0.0 * phi[0] - 0.5 * phi[1] + 0.5 * phi[2];
        wR[k][1] = 1.0 * phi[0] + 0.0 * phi[1] - 1.0 * phi[2];
        wR[k][2] = 0.0 * phi[0] + 0.5 * phi[1] + 0.5 * phi[2];
    }

    wR[0][0] = 1.0 * phi[0] - 1.5 * phi[1] + 0.5 * phi[2];
    wR[0][1] = 0.0 * phi[0] + 2.0 * phi[1] - 1.0 * phi[2];
    wR[0][2] = 0.0 * phi[0] - 0.5 * phi[1] + 0.5 * phi[2];


    // IL Sweep
    for( k=1; k<=M2; k++ )  {
        IL      =   d[k-1]*IL + 
                        wL[k-1][0] * u[0] + 
                        wL[k-1][1] * u[1] + 
                        wL[k-1][2] * u[2];      // wL(k,:)*u(iL);
        I[k]    =   I[k] + IL;
    }

    for( k=M2+1; k<=N-M2; k++ )  {
        IL      =   d[k-1]*IL + 
                        wL[k-1][0] * u[k-1] + 
                        wL[k-1][1] * u[k] + 
                        wL[k-1][2] * u[k+1];    // wL(k,:)*u(k+1+iC)
        I[k]    =   I[k] + IL;
    }

    for( k=N+1-M2; k<=N; k++ )  {
        IL      =   d[k-1]*IL + 
                        wL[k-1][0] * u[N-2] + 
                        wL[k-1][1] * u[N-1] + 
                        wL[k-1][2] * u[N];      // wL(k,:)*u(iR);
        I[k]    =   I[k] + IL;
    }




    // IR Sweep
    for( m=N; m>=N+1-M2; m-- )  {
        IR      =   d[m-1]*IR + 
                        wR[m-1][0] * u[N-2] + 
                        wR[m-1][1] * u[N-1] + 
                        wR[m-1][2] * u[N];      // wR(m,:)*u(iR);
        I[m-1]    =   I[m-1] + IR;
    }

    for( m=N-M2; m>=1+M2; m-- )  {
        IR      =   d[m-1]*IR + 
                        wR[m-1][0] * u[m-2] + 
                        wR[m-1][1] * u[m-1] + 
                        wR[m-1][2] * u[m];    // wR(m,:)*u(m+iC)
        I[m-1]    =   I[m-1] + IR;
    }

    for( m=M2; m>=1; m-- )  {
        IR      =   d[m-1]*IR + 
                        wR[m-1][0] * u[0] + 
                        wR[m-1][1] * u[1] + 
                        wR[m-1][2] * u[2];      // wR(m,:)*u(iL);
        I[m-1]    =   I[m-1] + IR;
    }



}



// Source term for the Maxwell solver
// TODO: Document this
void getSourceTerm( 
                    const dTensorBC2& qvals, 
                    dTensorBC2& auxvals )
{

    // Parameters
    double ion_mass = global_ini_params.get_ion_mass();
    double elc_mass = global_ini_params.get_elc_mass();

    const double cs_light_squared = global_ini_params.get_cs_light_squared();
    const double one_over_epsilon = cs_light_squared;

   
    // Loop over each point
    int numpts = auxvals.getsize(1);
    for (int i=1; i<=numpts; i++)
    {


        // Variables
        const double& M1_i     = qvals.get(i,_M1_i );
        const double& M2_i     = qvals.get(i,_M2_i );
        const double& M3_i     = qvals.get(i,_M3_i );
        
        const double& M1_e     = qvals.get(i,_M1_e );
        const double& M2_e     = qvals.get(i,_M2_e );
        const double& M3_e     = qvals.get(i,_M3_e );
        
        const double J1_i =  M1_i/ion_mass;
        const double J2_i =  M2_i/ion_mass;
        const double J3_i =  M3_i/ion_mass;
        const double J1_e = -M1_e/elc_mass;
        const double J2_e = -M2_e/elc_mass;
        const double J3_e = -M3_e/elc_mass;
        const double J1   = J1_i + J1_e;
        const double J2   = J2_i + J2_e;
        const double J3   = J3_i + J3_e;
        
        // Store the source terms in aux. storage for use by the
        // implicit Maxwell solver
        auxvals.set(i, 15, -J1*one_over_epsilon);
        auxvals.set(i, 16, -J2*one_over_epsilon);
        auxvals.set(i, 17, -J3*one_over_epsilon);
    }
}


// Function to impose (outflow) boundary conditions
// TODO: Other boundary conditions
// See [1] for details. This implements second order accurate outflow
// boundary conditions
void applyBdryCdn( 
        const int n, 
        const double I[], 
        const double u[], const double u_prev[], 
        const double alpha,     const double beta, 
        double& A_n, double& B_n, 
        const double xa, const double xb )
{


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
                Gamma_0*I[0] +      // Gamma_0 * I[0]
                Gamma_1*u[0] +      // Gamma_1 * u[0]
                Gamma_2*u_prev[0];  // Gamma_2 * u_prev[0]

    w_b     =   exp(-beta)*B_n + 
                Gamma_0*I[n] +      // Gamma_0 * I[n]     
                Gamma_1*u[n] +      // Gamma_1 * u[n]     
                Gamma_2*u_prev[n];  // Gamma_2 * u_prev[n]
                                        
    // This is where we solve the system
    A_n     =       ( (1.0-Gamma_0)*w_a + mu*Gamma_0*w_b ) / 
                    ( pow(1.0-Gamma_0,2.0) - pow(mu*Gamma_0,2.0) );
    B_n     =       ( (1.0-Gamma_0)*w_b + mu*Gamma_0*w_a ) / 
                    ( pow(1.0-Gamma_0,2.0) - pow(mu*Gamma_0,2.0) );                                


}





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
    double xa   = aux.get(1,13);         // Domain boundary
    double xb   = aux.get(n+1,13);
    double c    = global_ini_params.get_cs_light(); // Speed of light
    double h    = global_ini_params.get_dx();       // Grid size


    // Check if we have the right dt
    // FINESS sometimes takes two (or more) steps for each dt
    // The extra steps are usually for small dt
    double dt_requested = global_ini_params.get_max_dt();

    if( abs(dt - dt_requested) <= 1.0e-10 ) {

        // Some method-specific parameters (see [1] for details)
        // TODO: Move these to IniParams since they are used in several
        // places
        double beta     = 2.0;
        double nu       = beta*h/(c*dt);
        double alpha    = beta/(c*dt);
    
        // Some temporary storage
        int     ia      = 0;    // Loop inex
        int     ib      = 0;    // Loop inex
        double  u_nxt   = 0.0;  // Used to store solution at t^{n+1}
        // For applying boundary conditions
        double  A_n     = 0.0;
        double  B_n     = 0.0;
    

        // Some more temporary storage
        // TODO: Use auxiliary storage instead
        double I[n+1], u[n+1], u_prev[n+1];
        for (ia=0; ia<=n; ia++) {
            I[ia] = 0.0; u[ia] = 0.0; u_prev[ia] = 0.0;
        }

        // We will implement out Maxwell solver here
        // The Maxwell solver requires the solution of wave equations for
        // the scalar and vector potential. 
        // First, solve a wave equation for the scalar potential \phi
        // Read in potential values
        for (ia = 0; ia <= n; ia++) {
            u[ia] = aux.get(ia+1, 1);
            u_prev[ia] = aux.get(ia+1, 2);
        }
        // Implicit wave solver implementation
        // Step 1: 
        // Apply Green's function of modified Helmholtz operator 
        // (this is implemented in [1] as a fast convolution)
        // TODO: Use auxiliary storage for efficiency
        fastConvolve(I, u, n, nu);


        // Step 2: 
        // Apply boundary conitions
        // Note: only outflow boundary conditions available currently
        // TODO: Use auxiliary storage for efficiency
        A_n = aux.get(1, 14);        // Boundary constants
        B_n = aux.get(2, 14);
        applyBdryCdn( n, I, u, u_prev, 
                alpha, beta, A_n, B_n, xa, xb );


        // Step 3:
        // Compute solution for \phi at time t^{n+1}
        for (ia=0; ia<=n; ia++)  {
        // Solution at t^{n+1}
        u_nxt = -( pow(beta,2.0) - 2.0 )*u[ia] - 
                    u_prev[ia] + 
                    ( pow(beta,2.0)/2.0 )*( I[ia] + 
                        A_n*exp( -alpha*(aux.get(ia+1,13)-xa) ) + 
                        B_n*exp( -alpha*(xb-aux.get(ia+1,13)) ) );
        // Update solution at t^{n-1}
        aux.set(ia+1, 2, aux.get(ia+1, 1) );
        aux.set(ia+1, 1, u_nxt);

        // Save boundary constants
        aux.set(1, 14, A_n); 
        aux.set(2, 14, B_n);
    }




        double srcContribution = 0.0;
        getSourceTerm(qvals, aux);        // Evaluate source term


        // Next, solve a wave equation for the vector potential A
        // Start with component A_x (in 1D, these are uncoupled 
        // wave equations)
        // Read in potential values
        for (ia = 0; ia <= n; ia++) {
            u[ia] = aux.get(ia+1, 3);
            u_prev[ia] = aux.get(ia+1, 4);
            I[ia] = 0.0;
        }
        // Implicit wave solver implementation
        // Step 1: 
        // Apply Green's function of modified Helmholtz operator 
        // (this is implemented in [1] as a fast convolution)
        // TODO: Use auxiliary storage for efficiency
        fastConvolve(I, u, n, nu);


        // Step 2: 
        // Apply boundary conitions
        // Note: only outflow boundary conditions available currently
        // TODO: Use auxiliary storage for efficiency
        A_n = aux.get(3, 14);        // Boundary constants
        B_n = aux.get(4, 14);
        applyBdryCdn( n, I, u, u_prev, 
                alpha, beta, A_n, B_n, xa, xb );


        // Step 3:
        // Compute solution for \phi at time t^{n+1}
        for (ia=0; ia<=n; ia++)  {

        // Source contribution
        srcContribution = 0.0;
        for (ib=1; ib<=n+1; ib++)
            srcContribution = srcContribution + 
                aux.get(ib, 15) * exp( -alpha*abs( 
                        aux.get(ia+1,13)-aux.get(ib,13) ) );

        // Solution at t^{n+1}
        u_nxt = -( pow(beta,2.0) - 2.0 )*u[ia] - 
                    u_prev[ia] + 
                    ( pow(beta,2.0)/2.0 )*( I[ia] + 
                        (c*dt/beta)*srcContribution + 
                        A_n*exp( -alpha*(aux.get(ia+1,13)-xa) ) + 
                        B_n*exp( -alpha*(xb-aux.get(ia+1,13)) ) );
        // Update solution at t^{n-1}
        aux.set(ia+1, 5, aux.get(ia+1, 4) );
        aux.set(ia+1, 4, aux.get(ia+1, 3) );
        aux.set(ia+1, 3, u_nxt);

        // Save boundary constants
        aux.set(3, 14, A_n); 
        aux.set(4, 14, B_n);
    }




        // Next, we have component A_y
        // Read in potential values
        for (ia = 0; ia <= n; ia++) {
            u[ia] = aux.get(ia+1, 6);
            u_prev[ia] = aux.get(ia+1, 7);
            I[ia] = 0.0;
        }
        // Implicit wave solver implementation
        // Step 1: 
        // Apply Green's function of modified Helmholtz operator 
        // (this is implemented in [1] as a fast convolution)
        // TODO: Use auxiliary storage for efficiency
        fastConvolve(I, u, n, nu);


        // Step 2: 
        // Apply boundary conitions
        // Note: only outflow boundary conditions available currently
        // TODO: Use auxiliary storage for efficiency
        A_n = aux.get(5, 14);        // Boundary constants
        B_n = aux.get(6, 14);
        applyBdryCdn( n, I, u, u_prev, 
                alpha, beta, A_n, B_n, xa, xb );


        // Step 3:
        // Compute solution for \phi at time t^{n+1}
        for (ia=0; ia<=n; ia++)  {

        // Source contribution
        srcContribution = 0.0;
        for (ib=1; ib<=n+1; ib++)
            srcContribution = srcContribution + 
                aux.get(ib, 16) * exp( -alpha*abs( 
                        aux.get(ia+1,13)-aux.get(ib,13) ) );


        // Solution at t^{n+1}
        u_nxt = -( pow(beta,2.0) - 2.0 )*u[ia] - 
                    u_prev[ia] + 
                    ( pow(beta,2.0)/2.0 )*( I[ia] + 
                        (c*dt/beta)*srcContribution + 
                        A_n*exp( -alpha*(aux.get(ia+1,13)-xa) ) + 
                        B_n*exp( -alpha*(xb-aux.get(ia+1,13)) ) );
        // Update solution at t^{n-1}
        aux.set(ia+1, 8, aux.get(ia+1, 7) );
        aux.set(ia+1, 7, aux.get(ia+1, 6) );
        aux.set(ia+1, 6, u_nxt);

        // Save boundary constants
        aux.set(5, 14, A_n); 
        aux.set(6, 14, B_n);
    }






        // Finally, we have component A_z
        // Read in potential values
        for (ia = 0; ia <= n; ia++) {
            u[ia] = aux.get(ia+1, 9);
            u_prev[ia] = aux.get(ia+1, 10);
            I[ia] = 0.0;
        }
        // Implicit wave solver implementation
        // Step 1: 
        // Apply Green's function of modified Helmholtz operator 
        // (this is implemented in [1] as a fast convolution)
        // TODO: Use auxiliary storage for efficiency
        fastConvolve(I, u, n, nu);


        // Step 2: 
        // Apply boundary conitions
        // Note: only outflow boundary conditions available currently
        // TODO: Use auxiliary storage for efficiency
        A_n = aux.get(7, 14);        // Boundary constants
        B_n = aux.get(8, 14);
        applyBdryCdn( n, I, u, u_prev, 
                alpha, beta, A_n, B_n, xa, xb );


        // Step 3:
        // Compute solution for \phi at time t^{n+1}
        for (ia=0; ia<=n; ia++)  {

        // Source contribution
        srcContribution = 0.0;
        for (ib=1; ib<=n+1; ib++)
            srcContribution = srcContribution + 
                aux.get(ib, 16) * exp( -alpha*abs( 
                        aux.get(ia+1,13)-aux.get(ib,13) ) );

        // Solution at t^{n+1}
        u_nxt = -( pow(beta,2.0) - 2.0 )*u[ia] - 
                    u_prev[ia] + 
                    ( pow(beta,2.0)/2.0 )*( I[ia] + 
                        (c*dt/beta)*srcContribution + 
                        A_n*exp( -alpha*(aux.get(ia+1,13)-xa) ) + 
                        B_n*exp( -alpha*(xb-aux.get(ia+1,13)) ) );
        // Update solution at t^{n-1}
        aux.set(ia+1, 11, aux.get(ia+1, 10) );
        aux.set(ia+1, 10, aux.get(ia+1, 9) );
        aux.set(ia+1, 9, u_nxt);

        // Save boundary constants
        aux.set(7, 14, A_n); 
        aux.set(8, 14, B_n);
    }




    // Now, we evaluate the fields from the potentials
    // We use second order (one-sided) approximation of the time derivative
    for (ia=1; ia<=n+1; ia++) {
    // We start with the electric fields
    // E_x = -partial_t A_x - grad_x \phi
        if( ia>1 && ia<n+1 )
            qvals.set(ia, _E1, 
                -(1.5*aux.get(ia, 3)/dt - 2.0*aux.get(ia, 4)/dt 
                    + 0.5*aux.get(ia, 5)/dt) - 
                    ( aux.get(ia+1, 1) - aux.get(ia-1, 1) )/(2.0*h) );
        else if( ia == 1 )
            qvals.set(ia, _E1, 
                -(1.5*aux.get(ia, 3)/dt - 2.0*aux.get(ia, 4)/dt 
                    + 0.5*aux.get(ia, 5)/dt) - 
                    ( -1.5*aux.get(ia, 1) + 2.0*aux.get(ia+1, 1) - 
                      0.5*aux.get(ia+2, 1) )/h );
        else if( ia == n+1 )
            qvals.set(ia, _E1, 
                -(1.5*aux.get(ia, 3)/dt - 2.0*aux.get(ia, 4)/dt 
                    + 0.5*aux.get(ia, 5)/dt) -
                    ( 1.5*aux.get(ia, 1) - 2.0*aux.get(ia-1, 1) + 
                      0.5*aux.get(ia-2, 1) )/h );
    // E_y = -partial_t A_y
        qvals.set(ia, _E2, 
                -(1.5*aux.get(ia, 6)/dt - 2.0*aux.get(ia, 7)/dt 
                    + 0.5*aux.get(ia, 8)/dt) );

    // E_z = -partial_t A_z
        qvals.set(ia, _E3, 
                -(1.5*aux.get(ia, 9)/dt - 2.0*aux.get(ia, 10)/dt 
                    + 0.5*aux.get(ia, 11)/dt) );


    // Now, we set the magnetic fields
    // For the 1D problem, the B_x field is unaltered

    // B_y = -partial_x A_z
    // B_z =  partial_x A_y
        if( ia>1 && ia<n+1 )    {
            qvals.set(ia, _B2, 
                -( aux.get(ia+1, 9) - aux.get(ia-1, 9) )/(2.0*h) );
            qvals.set(ia, _B3, 
                 ( aux.get(ia+1, 6) - aux.get(ia-1, 6) )/(2.0*h) );
        }
        else if( ia == 1 )      {
            qvals.set(ia, _B2, 
                   -( -1.5*aux.get(ia, 9) + 2.0*aux.get(ia+1, 9) - 
                      0.5*aux.get(ia+2, 9) )/h );
            qvals.set(ia, _B3, 
                    ( -1.5*aux.get(ia, 6) + 2.0*aux.get(ia+1, 6) - 
                      0.5*aux.get(ia+2, 6) )/h );
        }
        else if( ia == n+1 )    {
            qvals.set(ia, _B2, 
                   -( 1.5*aux.get(ia, 9) - 2.0*aux.get(ia-1, 9) + 
                      0.5*aux.get(ia-2, 9) )/h );
            qvals.set(ia, _B3, 
                    ( 1.5*aux.get(ia, 6) - 2.0*aux.get(ia-1, 6) + 
                      0.5*aux.get(ia-2, 6) )/h );
        }


    }
    

    // TODO: Does Dave need average values of the fields?


    }


}
