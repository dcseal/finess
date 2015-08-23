#include "tensors.h"
#include <iostream>
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
//
//
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

void AuxFunc(const dTensor1& xpts, dTensor2& auxvals)
{

    // The "do nothing" aux func

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {

        // Spatial location.
        double x = xpts.get(i);

        // Note: We have to also supply initial conditions 
        // here for the auxiliary fields
        // TODO: Include documentation on how this was set
        // We use  
        // \phi = 0
        // \phi_t = 0
        // A such that curl(A) = B and div(A) = 0
        // A_t = 0

        // For the implicit Maxwell solver, hold three time levels
        // for the scalar and vector potential
        // for the scalar potential \phi^{n+1}, \phi^n and 
        // \phi^{n-1}
        auxvals.set(i, 1, 0.0);   // \phi^n
        auxvals.set(i, 2, 0.0);   // \phi^{n-1}

        // For the vector potential, A_x
        auxvals.set(i, 3, 0.0);   // A_x^n
        auxvals.set(i, 4, 0.0);   // A_x^{n-1}
        auxvals.set(i, 5, 0.0);   // A_x^{n-2}

        // For the vector potential, A_y
        auxvals.set(i, 6, 0.0);   // A_y^n      
        auxvals.set(i, 7, 0.0);   // A_y^{n-1}
        auxvals.set(i, 8, 0.0);   // A_y^{n-2}

        // For the vector potential, A_z
        if( x < 6.283185307179586 ) {
            // Left half of domain
            auxvals.set(i, 9, -x+2.0*M_PI);   // A_z^n    
            auxvals.set(i, 10, -x+2.0*M_PI);  // A_z^{n-1}
            auxvals.set(i, 11, -x+2.0*M_PI);  // A_z^{n-2}
        }
        else    {
            // Right half of domain
            auxvals.set(i, 9, x-2.0*M_PI);   // A_z^n    
            auxvals.set(i, 10, x-2.0*M_PI);  // A_z^{n-1}
            auxvals.set(i, 11, x-2.0*M_PI);  // A_z^{n-2}
        }


        // We will also create an auxiliary work array for evaluating 
        // successive convoutions efficiently
        // This holds the particular solution, I
        // (See Sec. 4 in [1] for details)
        auxvals.set(i, 12, 0.0);


        // Also store the grid points - we will use this 
        // in applying boundary conditions
        auxvals.set(i, 13, x);

        // Store boundary constants required for (outflow) boundary 
        // conditions
        // TODO: Find a better solution for this
        auxvals.set(i, 14, 0.0);

        // Store Maxwell source terms obtained from the fluid 
        // solve here. The Maxwell solver will make use of this
        auxvals.set(i, 15, 0.0);
        auxvals.set(i, 16, 0.0);
        auxvals.set(i, 17, 0.0);
        auxvals.set(i, 18, 0.0);

    }

}
