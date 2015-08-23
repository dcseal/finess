#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "Components.h"     // Used for indexing into solution
#include "IniParams.h"
#include "gas05.h"          // for common 5 moment operations

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

    const int numpts    = xpts.getsize();
    const int meqn      = qvals.getsize(2);
    const double& gamma = global_ini_params.get_gamma();

    // Brio-Wu parameters
    //
    // See Eqns. (45)-(46) in
    // "Approximate Riemann solver for the two-fluid plasma model"
    // for a description of what these parameters are.
    // Note that we have a minor rescaling of the problem here!
    const double& mass_ratio     = global_ini_params.get_mass_ratio  ();
    const double& debye_length   = global_ini_params.get_debye_length();
    const double& larmor_radius  = global_ini_params.get_larmor_radius();

    // Derived parameters
    const double ion_mass       = 1.0;
    const double elc_mass       = ion_mass/mass_ratio;


    // left values
    //
    double const& left_rho_i   = 1.0*ion_mass;
    double const& left_u1_i    = 0.0;
    double const& left_u2_i    = 0.0;
    double const& left_u3_i    = 0.0;
    double const& left_press_i = 0.5;

    double const& left_rho_e   = 1.0*elc_mass;
    double const& left_u1_e    = 0.0;
    double const& left_u2_e    = 0.0;
    double const& left_u3_e    = 0.0;
    double const& left_press_e = 0.5;

    double const& left_B1      = 0.75;
    double const& left_B2      = 1.0;
    double const& left_B3      = 0.0;
    double const& left_E1      = 0.0;
    double const& left_E2      = 0.0;
    double const& left_E3      = 0.0;	    

    // right values
    //
    double const& rght_rho_i   = 0.125*ion_mass;
    double const& rght_u1_i    = 0.0;
    double const& rght_u2_i    = 0.0;
    double const& rght_u3_i    = 0.0;
    double const& rght_press_i = 0.05;

    double const& rght_rho_e   = 0.125*elc_mass;
    double const& rght_u1_e    = 0.0;
    double const& rght_u2_e    = 0.0;
    double const& rght_u3_e    = 0.0;
    double const& rght_press_e = 0.05;

    double const& rght_B1      = 0.75;
    double const& rght_B2      =-1.0;
    double const& rght_B3      = 0.0;
    double const& rght_E1      = 0.0;
    double const& rght_E2      = 0.0;
    double const& rght_E3      = 0.0;

    //
    // derived values
    //
    double const left_M1_i = left_rho_i*left_u1_i;
    double const left_M2_i = left_rho_i*left_u2_i;
    double const left_M3_i = left_rho_i*left_u3_i;
    //
    double const left_M1_e = left_rho_i*left_u1_e;
    double const left_M2_e = left_rho_i*left_u2_e;
    double const left_M3_e = left_rho_i*left_u3_e;
    //
    double const left_energy_i = left_press_i/(gamma-1.0e0) + 0.5e0*left_rho_i
        *(left_u1_i*left_u1_i + left_u2_i*left_u2_i + left_u3_i*left_u3_i);
    double const left_energy_e = left_press_e/(gamma-1.0e0) + 0.5e0*left_rho_e
        *(left_u1_e*left_u1_e + left_u2_e*left_u2_e + left_u3_e*left_u3_e);
    const double  left_entropy_i = get_entropy_per_vol05(left_rho_i,left_press_i);
    const double  left_entropy_e = get_entropy_per_vol05(left_rho_e,left_press_e);
    //
    double const rght_M1_i = rght_rho_i*rght_u1_i;
    double const rght_M2_i = rght_rho_i*rght_u2_i;
    double const rght_M3_i = rght_rho_i*rght_u3_i;
    //
    double const rght_M1_e = rght_rho_i*rght_u1_e;
    double const rght_M2_e = rght_rho_i*rght_u2_e;
    double const rght_M3_e = rght_rho_i*rght_u3_e;
    //
    double const rght_energy_i = rght_press_i/(gamma-1.0e0) + 0.5e0*rght_rho_i
        *(rght_u1_i*rght_u1_i + rght_u2_i*rght_u2_i + rght_u3_i*rght_u3_i);
    double const rght_energy_e = rght_press_e/(gamma-1.0e0) + 0.5e0*rght_rho_e
        *(rght_u1_e*rght_u1_e + rght_u2_e*rght_u2_e + rght_u3_e*rght_u3_e);
    const double  rght_entropy_i = get_entropy_per_vol05(rght_rho_i,rght_press_i);
    const double  rght_entropy_e = get_entropy_per_vol05(rght_rho_e,rght_press_e);

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i);

        if( x < 6.283185307179586 )
        {
            qvals.set(i,_rho_i,  left_rho_i );
            qvals.set(i,_M1_i ,  left_M1_i );
            qvals.set(i,_M2_i ,  left_M2_i );
            qvals.set(i,_M3_i ,  left_M3_i );
            qvals.set(i,_N_i  ,  left_energy_i );

            qvals.set(i,_rho_e,  left_rho_e );
            qvals.set(i,_M1_e ,  left_M1_e );
            qvals.set(i,_M2_e ,  left_M2_e );
            qvals.set(i,_M3_e ,  left_M3_e );
            qvals.set(i,_N_e  ,  left_energy_e );

            qvals.set(i,_B1   , left_B1 );
            qvals.set(i,_B2   , left_B2 );
            qvals.set(i,_B3   , left_B3 );
            qvals.set(i,_E1   , left_E1 );
            qvals.set(i,_E2   , left_E2 );
            qvals.set(i,_E3   , left_E3 );

            // B-field, (psi) and then E-field (phi) cleaning
//          if(meqn<_psi) continue;
//          qvals.set(i,_psi  , 0. );
//          if(meqn<_phi) continue;
//          qvals.set(i,_phi  , 0. );

//          // entropy tracking
//          if(meqn >=_entropy_i)
//          {
//              qvals.set(i,_entropy_i, left_entropy_i );
//              qvals.set(i,_entropy_e, left_entropy_e );
//          }
        }
        else
        {
            qvals.set(i,_rho_i,  rght_rho_i );
            qvals.set(i,_M1_i ,  rght_M1_i );
            qvals.set(i,_M2_i ,  rght_M2_i );
            qvals.set(i,_M3_i ,  rght_M3_i );
            qvals.set(i,_N_i  ,  rght_energy_i );

            qvals.set(i,_rho_e,  rght_rho_e );
            qvals.set(i,_M1_e ,  rght_M1_e );
            qvals.set(i,_M2_e ,  rght_M2_e );
            qvals.set(i,_M3_e ,  rght_M3_e );
            qvals.set(i,_N_e  , rght_energy_e );

            qvals.set(i,_B1   , rght_B1 );
            qvals.set(i,_B2   , rght_B2 );
            qvals.set(i,_B3   , rght_B3 );
            qvals.set(i,_E1   , rght_E1 );
            qvals.set(i,_E2   , rght_E2 );
            qvals.set(i,_E3   , rght_E3 );

//          if(meqn<_psi) continue;
//          qvals.set(i,_psi  , 0. );
//          if(meqn<_phi) continue;
//          qvals.set(i,_phi  , 0. );

//          if(meqn >=_entropy_i)
//          {
//              qvals.set(i,_entropy_i, rght_entropy_i );
//              qvals.set(i,_entropy_e, rght_entropy_e );
//          }
        }
    }
}
