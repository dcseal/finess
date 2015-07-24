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
        dTensorBC2& aux, 
        const double nu   )  
{

    // Basic parameters
    int         N   = aux.getsize(1) - 1;   // In [1], the domain is 
                                            // discretized using n+1 points
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

    // Store solution u^n here
    double u[N+1];
    for( k=0; k<=N; k++ )
        u[k] = aux.get(k+1, 1);       // Copy from aux. var

    // Initialize integral solution, I
    double I[N+1];
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




    // Copy computed I back to aux. variable structure
    for( k=0; k<N; k++ )
        aux.set( k+1, 3, I[k] );


}


