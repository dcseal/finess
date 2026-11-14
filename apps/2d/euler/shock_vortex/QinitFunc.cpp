#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
//      2D Euler equations.
//
//      Q(1) = density
//      Q(2) = rho u^1
//      Q(3) = rho u^2
//      Q(4) = rho u^3
//      Q(5) = Energy
//
// Input:
//
//    xpts( 1:numpts, 1:2 )      - The x, and y-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: $FINESS/lib/2d/blanks/QinitFunc
// Free-stream (upstream) state, shared with SetBndValues.cpp.  The upstream
// gas has density rho0 and pressure p0 (default rho0 = p0 = 1.0, as in the
// paper), and it moves to the right with Mach number MS relative to its
// sound speed c0 = sqrt(gamma*p0/rho0).  The shock is stationary at x = x0.
void ShockVortexUpstream(double& rho, double& u1, double& press)
{
    const double gamma = global_ini_params.get_gamma();
    const double MS    = 1.5;
    const double rho0  = 1.0;
    const double p0    = 1.0;
    rho   = rho0;
    press = p0;
    u1    = MS*sqrt(gamma*p0/rho0);
}

void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

// See 5.7. Shock-vortex interaction in "A posteriori subcell limiting of the
// discontinuous Galerkin finite element method for hyperbolic conservation
// laws," Dumbser et al, JCP 2014.
//
// A stationary normal shock sits at x = x0.  A vortex, centered at (xc,yc)
// in the upstream flow, is carried into the shock.  The vortex has angular
// velocity
//
//      v_phi(r) = v_m r/a                          for     r <= a,
//      v_phi(r) = v_m a/(a^2-b^2) * (r - b^2/r)    for a < r <= b,
//      v_phi(r) = 0                                for b < r,
//
// with v_m = M_v*c0, where c0 = sqrt(gamma*p0/rho0) is the upstream sound
// speed.  Its temperature T = p/rho (gas constant R = 1) solves
// dT/dr = (gamma-1)/gamma * v_phi^2/r with T(b) = T0, and density and
// pressure follow from T isentropically.

    // Gas consant
    const double gamma = global_ini_params.get_gamma();
    const double x0    = global_ini_params.get_x0();

    // Parameters for this probelm
    const double a    = 0.075; 
    const double b    = 0.175; 
    const double M_v  = 0.7;
    const double xc   = 0.25;
    const double yc   = 0.5;

    // Upstream state, and its Mach number MS
    double rho0, u0, p0;
    ShockVortexUpstream(rho0, u0, p0);
    const double T0  = p0/rho0;
    const double c0  = sqrt(gamma*p0/rho0);
    const double MS  = u0/c0;
    const double v_m = M_v*c0;

    // Downstream state (Rankine-Hugoniot conditions for a stationary shock)
    const double M2   = MS*MS;
    const double rhoR = rho0*(gamma+1.0)*M2/((gamma-1.0)*M2 + 2.0);
    const double pR   = p0*(1.0 + 2.0*gamma/(gamma+1.0)*(M2-1.0));
    const double uR   = u0*rho0/rhoR;

    // Integral of v_phi^2/s from r to b (exact; needs a < r <= b)
    const double K = v_m*a/(a*a - b*b);
    auto outer_integral = [&](double r)
    {
        return K*K*( 0.5*(b*b - r*r) - 2.0*b*b*log(b/r)
                   + 0.5*b*b*b*b*(1.0/(r*r) - 1.0/(b*b)) );
    };
    const double Ia = outer_integral(a);

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        double rho,u1,u2,u3,press;
        u3 = 0.0;

        if( x >= x0 )
        {
            rho   = rhoR;
            u1    = uR;
            u2    = 0.0;
            press = pR;
        }
        else
        {
            // Distance to center
            const double r = sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) );

            // Angular velocity and the integral of v_phi^2/s from r to b
            double v_phi, I;
            if( r <= a )
            {
                v_phi = v_m * r / a;
                I     = Ia + 0.5*v_m*v_m*(a*a - r*r)/(a*a);
            }
            else if( r <= b )
            {
                v_phi = K * (r - b*b/r);
                I     = outer_integral(r);
            } 
            else
            {
                v_phi = 0.;
                I     = 0.;
            }

            const double T = T0 - (gamma-1.0)/gamma*I;
            rho   = rho0*pow(T/T0, 1.0/(gamma-1.0));
            press = p0  *pow(T/T0, gamma/(gamma-1.0));

            // Counter-clockwise rotation about (xc,yc), added to the upstream
            // flow: (u1,u2) = (u0,0) + v_phi*(-sin(theta), cos(theta)), where
            // theta is the polar angle of (x,y) about the vortex center.
            const double sin_theta = (r > 0.0) ? (y-yc)/r : 0.0;
            const double cos_theta = (r > 0.0) ? (x-xc)/r : 0.0;
            u1 = u0 - v_phi*sin_theta;
            u2 =      v_phi*cos_theta;
        }

        // Comptue derived quantities and store each of the variables
        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
