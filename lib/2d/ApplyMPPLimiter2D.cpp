#include <cmath>
#include<iostream>
#include "stdio.h"
#include "dog_math.h"
#include "IniParams.h"
#include "constants.h"
#include "tensors.h"
using namespace std;

// MPP limiter.  This needs to be written individually for each application
// that attempts to use it.
//
// This limiter works by considering a high order flux, (fHat, gHat), and
// low-order flux, (fLF, gLF), and then creates a new flux of the form:
//
//     fhat = theta fhat + (1-theta) fLF
//     ghat = theta ghat + (1-theta) gLg
//
// so that the conservative update, 
//
// q = q - dt/dx( f_{i+1/2} - f_{i-1/2} ) - dt/dy( g_{j+1/2} - g_{j-1/2} )
//
// is conservative.
//
// See: "http://arxiv.org/abs/1411.0328" and references therein for more
// details.
void ApplyMPPLimiter2D( 
        const double dt, const dTensorBC3& q, 
        const dTensorBC3& fLF, const dTensorBC3& gLF,
        dTensorBC3& fHat, dTensorBC3& gHat )
{
    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();

    printf("WARNING: no positivity limiter has been implemented for this application\n");

}
