#include <cmath>
#include <iostream>
#include "dogdefs.h"
#include "stdio.h"
#include "IniParams.h"
using namespace std;

// MPP limiter.  This needs to be written individually for each application
// that attempts to use it.
//
// This limiter works by considering a high order flux, fHat, and
// low-order flux, fLF, and then creates a new flux of the form:
//
//     fhat = theta fhat + (1-theta) fLF
//
// so that the conservative update, 
//
// q = q - dt/dx( f_{i+1/2} - f_{i-1/2} )
//
// retains positivity of the density and pressure.
//
// See: "http://arxiv.org/abs/1411.0328" and references therein for more
// details.
void ApplyMPPLimiter1D( const double dt, const dTensorBC2& q, const dTensorBC2& fLF, dTensorBC2& fHat )
{
    // Parameters for the current grid
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);

    printf("WARNING: no positivity limiter has been implemented for this application\n");

}
