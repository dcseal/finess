#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-required routine that defines the initial conditions for the
// problem.
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{

    // Read in some problem parameters
    const int numpts    = xpts.getsize();

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i);

        // For now, we will hard-code the initial condition
        // We will simulate an exponential (two-way) wave problem
        // u(x,t) = exp( -25 * (x-1/4+ct)^2 ) + exp( -25 * (x+1/4-ct)^2 )
        // Hence, initial condition is 
        // u(x,0) =  exp( -25 * (x-1/4)^2 ) + exp( -25 * (x+1/4)^2 )
        // Note: Since this is a second order equation, we have two initial 
        // conditions. However, in an effort to not make too many changes to 
        // the structure of FINESS, we will hard-code this in the 
        // auxiliary variables 
        // // (see AuxFunc.cpp in ../lib/)
        qvals.set(i, 1,     exp( -25.0*pow(x+0.25, 2.0) ) + 
                exp( -25.0*pow(x-0.25, 2.0) )       );

    }
}
