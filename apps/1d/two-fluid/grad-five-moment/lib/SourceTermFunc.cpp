#include "tensors.h"
#include "Maxwell.h"
#include "Components.h"
#include "IniParams.h"

// Source term function psi in the hyperbolic balance law:
//
//            q_t + f_x = psi,
//
// Input:
//
//       xpts( 1:numpts )           - The x-coordinates for a list of points
//      qvals( 1:numpts, 1:meqn )   - Solution q at each point.
//    auxvals( 1:numpts, 1:maux )   - The auxilary function at each point.
//
// Output:
//
//     psi( 1:numpts, 1:meqn )      - The source term evaluated at each point.
//
// See also: FluxFunc.
void SourceTermFunc(const dTensor1& xpts, 
                    const dTensor2& qvals, 
                    const dTensor2& auxvals,
                    dTensor2& source )
{

    // Parameters
    //
    // We use the nondimensionalization found in "Entropy Stable Numerical
    // Schemes for Two-Fluid Plasma Equation", J. Sci. Comput. (2012)
    //
    // Note that this nondimensionalization implies that rho_i = n_i (ion
    // density is same as ion number density) because the ion mass is assumed
    // to be unity.

    const double mass_ratio     = global_ini_params.get_mass_ratio  ();
    const double debye_length   = global_ini_params.get_debye_length();
    const double larmor_radius  = global_ini_params.get_larmor_radius();

    //const double mass_ratio     = global_ini_params.get_mass_ratio  ();
    const double debye = global_ini_params.get_debye_length();
    //const double larmor_radius  = global_ini_params.get_larmor_radius();


    // Derived parameters
    const double ion_mass       = 1.0;
    const double elc_mass       = ion_mass/mass_ratio;

    // Maxwell parameters
    const double cs_light_squared = global_ini_params.get_cs_light_squared();
    const double cp_speed_squared = cs_light_squared*global_ini_params.get_cc_sqd();
    const double one_over_epsilon = 1.0 / (debye_length*debye_length*larmor_radius);

    // Loop over each point
    int numpts = xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i);

        // Variables
//      const double& rho_i    = qvals.get(i,_rho_i);
//      const double& M1_i     = qvals.get(i,_M1_i );
//      const double& M2_i     = qvals.get(i,_M2_i );
//      const double& M3_i     = qvals.get(i,_M3_i );
//      
//      const double& rho_e    = qvals.get(i,_rho_e);
//      const double& M1_e     = qvals.get(i,_M1_e );
//      const double& M2_e     = qvals.get(i,_M2_e );
//      const double& M3_e     = qvals.get(i,_M3_e );
//      
//      const double& B1       = qvals.get(i,_B1   );
//      const double& B2       = qvals.get(i,_B2   );
//      const double& B3       = qvals.get(i,_B3   );
//      const double& E1       = qvals.get(i,_E1   );
//      const double& E2       = qvals.get(i,_E2   );
//      const double& E3       = qvals.get(i,_E3   );

//      const double ion_charge_dens =  rho_i;
//      const double J1_i =  M1_i/ion_mass;
//      const double J2_i =  M2_i/ion_mass;
//      const double J3_i =  M3_i/ion_mass;
//      const double J1_e = -M1_e/elc_mass;
//      const double J2_e = -M2_e/elc_mass;
//      const double J3_e = -M3_e/elc_mass;
//      const double J1   = J1_i + J1_e;
//      const double J2   = J2_i + J2_e;
//      const double J3   = J3_i + J3_e;
//      
//      // Source term
//      source.set(i,_rho_i,   0.0 );
//      source.set(i,_M1_i ,   ( rho_i*E1 + (J2_i*B3 - J3_i*B2) ) / larmor_radius );
//      source.set(i,_M2_i ,   ( rho_i*E2 + (J3_i*B1 - J1_i*B3) ) / larmor_radius );
//      source.set(i,_M3_i ,   ( rho_i*E3 + (J1_i*B2 - J2_i*B1) ) / larmor_radius );
//      source.set(i,_N_i  ,   (  J1_i*E1 +  J2_i*E2 + J3_i*E3  ) / larmor_radius );

//      source.set(i,_rho_e,   0.0 );
//      source.set(i,_M1_e ,   -rho_e*(E1 + (J2_e*B3 - J3_e*B2) )*mass_ratio / larmor_radius );
//      source.set(i,_M2_e ,   -rho_e*(E2 + (J3_e*B1 - J1_e*B3) )*mass_ratio / larmor_radius );
//      source.set(i,_M3_e ,   -rho_e*(E3 + (J1_e*B2 - J2_e*B1) )*mass_ratio / larmor_radius );
//      source.set(i,_N_e  ,   -rho_e*(J1_e*E1 + J2_e*E2 + J3_e*E3)*mass_ratio / larmor_radius );

//      source.set(i,_B1   ,   0.0 );
//      source.set(i,_B2   ,   0.0 );
//      source.set(i,_B3   ,   0.0 );
//      source.set(i,_E1   ,   -J1*one_over_epsilon);
//      source.set(i,_E2   ,   -J2*one_over_epsilon);
//      source.set(i,_E3   ,   -J3*one_over_epsilon);

        // ---- This (below here) was copied directly from DogPack ---- //

        // Variables
        const double rho_i    = qvals.get(i,1);
        const double u1_i     = qvals.get(i,2)/rho_i;
        const double u2_i     = qvals.get(i,3)/rho_i;
        const double u3_i     = qvals.get(i,4)/rho_i;

        const double rho_e    = qvals.get(i,6);
        const double u1_e     = qvals.get(i,7)/rho_e;
        const double u2_e     = qvals.get(i,8)/rho_e;
        const double u3_e     = qvals.get(i,9)/rho_e;

        const double B1       = qvals.get(i,11);
        const double B2       = qvals.get(i,12);
        const double B3       = qvals.get(i,13);
        const double E1       = qvals.get(i,14);
        const double E2       = qvals.get(i,15);
        const double E3       = qvals.get(i,16);

        // Source term
        source.set(i,1,   0.0 );
        source.set(i,2,   rho_i*(E1 + u2_i*B3 - u3_i*B2)/larmor_radius );
        source.set(i,3,   rho_i*(E2 + u3_i*B1 - u1_i*B3)/larmor_radius );
        source.set(i,4,   rho_i*(E3 + u1_i*B2 - u2_i*B1)/larmor_radius );
        source.set(i,5,   rho_i*(u1_i*E1 + u2_i*E2 + u3_i*E3)/larmor_radius );
        source.set(i,6,   0.0 );
        source.set(i,7,  -rho_e*(E1 + u2_e*B3 - u3_e*B2)*mass_ratio/larmor_radius );
        source.set(i,8,  -rho_e*(E2 + u3_e*B1 - u1_e*B3)*mass_ratio/larmor_radius );
        source.set(i,9,  -rho_e*(E3 + u1_e*B2 - u2_e*B1)*mass_ratio/larmor_radius );
        source.set(i,10, -rho_e*(u1_e*E1 + u2_e*E2 + u3_e*E3)*mass_ratio/larmor_radius );
        source.set(i,11,  0.0 );
        source.set(i,12,  0.0 );
        source.set(i,13,  0.0 );
        source.set(i,14,  (mass_ratio*rho_e*u1_e-rho_i*u1_i)/(debye*debye*larmor_radius) );
        source.set(i,15,  (mass_ratio*rho_e*u2_e-rho_i*u2_i)/(debye*debye*larmor_radius) );
        source.set(i,16,  (mass_ratio*rho_e*u3_e-rho_i*u3_i)/(debye*debye*larmor_radius) );

/*
 * TODO - reinclude this section! -DS 7/9/2015

        // Check for hyperbolic cleaning corrections
        const int& meqn = source.getsize(2);
        if(_psi > meqn) continue;
        source.set(i, _psi,  0. );

        if(_phi > meqn) continue;
        if( charge_density_coef )
        {
          const double phi = qvals.get(i,_phi);
          source.set(i,_phi, charge_density_coef*charge_density - eps_Ecp_func(phi,time) );
        }
        else
        { source.set(i,_phi, 0.0); }

*/

// Extra entropy equations?
//      if(_entropy_i > meqn) continue;
//      source.set(i,_entropy_i, 0.0);
//      source.set(i,_entropy_e, 0.0);

    }
}
