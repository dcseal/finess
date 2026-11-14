#include <cmath>
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"

// This is a user-supplied routine that sets the boundary conditions for the
// shock-vortex interaction problem:
//
//      left:         supersonic inflow (the upstream state of QinitFunc.cpp),
//      right:        zeroth-order extrapolation (outflow),
//      top, bottom:  reflective (solid) walls.
//
void SetBndValues( StateVars& Q )
{

    dTensorBC3& q    = Q.ref_q();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();

    // Gas constant
    const double gamma = global_ini_params.get_gamma();

    // ********************************************************************* //
    // LEFT BOUNDARY (supersonic inflow)
    // ********************************************************************* //
    void ShockVortexUpstream(double& rho, double& u1, double& press);
    double rho, u1, press;
    ShockVortexUpstream(rho, u1, press);
    const double energy = press/(gamma-1.0e0) + 0.5e0*rho*u1*u1;
    for (int i=1-mbc; i<=0; i++)
    for (int j=1; j<=my; j++)
    {
        q.set(i,j,1, rho    );
        q.set(i,j,2, rho*u1 );
        q.set(i,j,3, 0.0    );
        q.set(i,j,4, 0.0    );
        q.set(i,j,5, energy );
    }
    // ********************************************************************* //

    // ********************************************************************* //
    // RIGHT BOUNDARY (Zeroth-order extrapolation)
    // ********************************************************************* //
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my; j++)
    for (int m=1; m<=meqn; m++)
    {
        q.set(i,j,m, q.get(mx,j,m) );
    }
    // ********************************************************************* //

    // ********************************************************************* //
    // BOTTOM AND TOP BOUNDARIES (reflective walls).  These run over the full
    // width, including the left/right ghost cells set above, so the corners
    // get filled in too.
    // ********************************************************************* //
    for (int i=1-mbc; i<=mx+mbc; i++)
    for (int k=1; k<=mbc; k++)
    {
        for (int m=1; m<=meqn; m++)
        {
            q.set(i, 1-k,  m, q.get(i, k,      m) );
            q.set(i, my+k, m, q.get(i, my+1-k, m) );
        }
        q.set(i, 1-k,  3, -q.get(i, 1-k,  3) );
        q.set(i, my+k, 3, -q.get(i, my+k, 3) );
    }
    // ********************************************************************* //

}

// Wrappers for main Euler library
void SetBndValuesX(StateVars& Q)
{ SetBndValues( Q ); }

void SetBndValuesY(StateVars& Q)
{ SetBndValues( Q ); }
