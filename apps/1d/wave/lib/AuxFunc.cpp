#include "tensors.h"
#include "IniParams.h"
#include <cmath>


using namespace std; 


// *REQUIRED*
//
// This is a user-required routine that defines the auxilary arrays when the
// parameter maux > 0.
//
// Each application is REQUIRED to define one of these.
//
//
// Input:
//
//    xpts( 1:numpts )             - The x-coordinates for a list of points
//
// Output:
//
//    auxvals( 1:numpts, 1:maux )  - The vector containg auxilary values.
//
// See also: QinitFunc.
void AuxFunc(const dTensor1& xpts, dTensor2& auxvals)
{

    // Read in some basic parameters
    const int numpts=xpts.getsize();
    const double dt = global_ini_params.get_initial_dt();
    const double c = global_ini_params.get_cs_light();


    for (int i=1; i<=numpts; i++)
    {

        // Spatial location.
        double x = xpts.get(i);

        // The implicit solver we use is a multi-step method
        // Current implementation is 2nd order in space and time
        // This requires us to store the solution at t^n and t^{n-1}
        //
        // For a description of the method used, see
        // [1]  Method of Lines Transpose: A Fast Implicit Wave Propagator
        //      Matthew F. Causley, Andrew J. Christlieb, Yaman Guclu, Eric Wolf
        //      arXiv:1306.6902
        //      2013
        //
        // Since the auxiliary variables are going to hold the solution, 
        // we will now provide initial conditions here
        // We will simulate an exponential (two-way) wave problem
        // u(x,t) = exp( -25 * (x-1/4+ct)^2 ) + exp( -25 * (x+1/4-ct)^2 )

        // Initialize u^n
        auxvals.set(i, 1,     exp( -25.0*pow(x+0.25, 2.0) ) + 
                exp( -25.0*pow(x-0.25, 2.0) )       );

        // Initialize u^{n-1}
        // To avoid taking a single step using an explicit method, we will 
        // set u(x, -dt)
        auxvals.set(i, 2,     exp( -25.0*pow(x+0.25+c*dt, 2.0) ) + 
                exp( -25.0*pow(x-0.25-c*dt, 2.0) )       );

        // We will also create an auxiliary work array for evaluating 
        // successive convoutions efficiently
        // This holds the particular solution, I
        // (See Sec. 4 in [1] for details)
        auxvals.set(i, 3, 0.0);

        // Also store the grid points - we will use this 
        // in applying boundary conditions
        auxvals.set(i, 4, x);
    }


    // For appliying outflow boundary conditions, we need a time 
    // history of two boundary constants
    // We will store these in the aux vars.
    // In particular, A_n = I[0] and B_n = I[n+2]
    auxvals.set(0, 3, 0.0);
    auxvals.set(numpts+1, 3, 0.0);


}
