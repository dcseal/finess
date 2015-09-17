#include "dogdefs.h"
#include "StateVars.h"
#include "IniParams.h"

// Function that is called after initial condition
//
// The purpose of this routine is to insert a single delta function into the
// lower-left corner of the domain.  The reason for doing this, is to make
// sure that the correct initial conditions are set.
//
// See also: Qinit.cpp
void AfterQinit( StateVars& Q )
{

    dTensorBC3& q   = Q.ref_q();
    dTensorBC3& aux = Q.ref_aux();

    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

  double energy = 0.244816/(dx*dy);
//    double energy = 0.979264/(dx*dy);

    q.set(1,1,5, energy );

}
