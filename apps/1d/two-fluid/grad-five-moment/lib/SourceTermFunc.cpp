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

//  2D version of source term function
//  void SourceTermFunc( const dTensor2& xpts, const dTensor2& qvals,
//       const dTensor2& auxvals, dTensor2& source)


    // Parameters
    double ion_mass = global_ini_params.get_ion_mass();
    double elc_mass = global_ini_params.get_elc_mass();

    const double cs_light_squared = global_ini_params.get_cs_light_squared();
    const double cp_speed_squared = cs_light_squared*global_ini_params.get_cc_sqd();
    const double one_over_epsilon = cs_light_squared;

    double charge_density_coef = cp_speed_squared*one_over_epsilon;
    if( !global_ini_params.get_clean_E_field() ) charge_density_coef = 0.;
   
    // TODO - fix this so the correct time gets passed in!  -DS
    // const double time = DogSolver::get_time_hack();
    
    // Loop over each point
    int numpts = xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i);
//      double y = xpts.get(i,2);

        // Variables
        const double& rho_i    = qvals.get(i,_rho_i);
        const double& M1_i     = qvals.get(i,_M1_i );
        const double& M2_i     = qvals.get(i,_M2_i );
        const double& M3_i     = qvals.get(i,_M3_i );
        
        const double& rho_e    = qvals.get(i,_rho_e);
        const double& M1_e     = qvals.get(i,_M1_e );
        const double& M2_e     = qvals.get(i,_M2_e );
        const double& M3_e     = qvals.get(i,_M3_e );
        
        const double& B1       = qvals.get(i,_B1   );
        const double& B2       = qvals.get(i,_B2   );
        const double& B3       = qvals.get(i,_B3   );
        const double& E1       = qvals.get(i,_E1   );
        const double& E2       = qvals.get(i,_E2   );
        const double& E3       = qvals.get(i,_E3   );

        const double ion_charge_dens =  rho_i/ion_mass;
        const double elc_charge_dens = -rho_e/elc_mass;
        const double charge_density  = ion_charge_dens + elc_charge_dens;
        const double J1_i =  M1_i/ion_mass;
        const double J2_i =  M2_i/ion_mass;
        const double J3_i =  M3_i/ion_mass;
        const double J1_e = -M1_e/elc_mass;
        const double J2_e = -M2_e/elc_mass;
        const double J3_e = -M3_e/elc_mass;
        const double J1   = J1_i + J1_e;
        const double J2   = J2_i + J2_e;
        const double J3   = J3_i + J3_e;
        
        // Source term
        source.set(i,_rho_i,   0.0 );
        source.set(i,_M1_i ,   ion_charge_dens*E1 + (J2_i*B3 - J3_i*B2));
        source.set(i,_M2_i ,   ion_charge_dens*E2 + (J3_i*B1 - J1_i*B3));
        source.set(i,_M3_i ,   ion_charge_dens*E3 + (J1_i*B2 - J2_i*B1));
        source.set(i,_N_i  ,   J1_i*E1 + J2_i*E2 + J3_i*E3);

        source.set(i,_rho_e,   0.0 );
        source.set(i,_M1_e ,   elc_charge_dens*E1 + (J2_e*B3 - J3_e*B2));
        source.set(i,_M2_e ,   elc_charge_dens*E2 + (J3_e*B1 - J1_e*B3));
        source.set(i,_M3_e ,   elc_charge_dens*E3 + (J1_e*B2 - J2_e*B1));
        source.set(i,_N_e  ,   (J1_e*E1 + J2_e*E2 + J3_e*E3));

        source.set(i,_B1   ,   0.0 );
        source.set(i,_B2   ,   0.0 );
        source.set(i,_B3   ,   0.0 );
        source.set(i,_E1   ,   -J1*one_over_epsilon);
        source.set(i,_E2   ,   -J2*one_over_epsilon);
        source.set(i,_E3   ,   -J3*one_over_epsilon);

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
