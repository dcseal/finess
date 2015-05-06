#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
//
// TODO - this routine needs some documentation! -DS
void SourceTermFunc( const dTensor2& xpts, const dTensor2& qvals,
     const dTensor2& auxvals, dTensor2& source)
{

//  if( global_ini_params.get_source_type_Emom() != SourceType::SIMPLE )
//  {
//      eprintf("You need to set source_term = 0 in parameters.ini");
//  }
    eprintf("You need to set a source term for this problem\n");

/*
    // Parameters
    double ion_mass = global_ini_params.get_ion_mass();
    double elc_mass = global_ini_params.get_elc_mass();
    const double one_over_epsilon = global_ini_params.get_cs_light_squared();
    const double cp_speed_squared = global_ini_params.get_cp_speed_squared();
    const double charge_density_coef = cp_speed_squared*one_over_epsilon;
    if(!maxwellParams.get_clean_E_field()) charge_density_coef = 0.;
   
    // TODO - FIGURE OUT HOW TO PASS TIME IN HERE! -DS
    // const double time = DogSolver::get_time_hack();
    const double time = 0.0;
    
    int numpts=xpts.getsize(1);
    // Loop over each point
    for (int i=1; i<=numpts; i++)
    {

        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

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

        const double ion_charge_dens = rho_i/ion_mass;
        const double elc_charge_dens = -rho_e/elc_mass;
        const double charge_density = ion_charge_dens + elc_charge_dens;
        const double J1_i = M1_i/ion_mass;
        const double J2_i = M2_i/ion_mass;
        const double J3_i = M3_i/ion_mass;
        const double J1_e = -M1_e/elc_mass;
        const double J2_e = -M2_e/elc_mass;
        const double J3_e = -M3_e/elc_mass;
        const double J1 = J1_i + J1_e;
        const double J2 = J2_i + J2_e;
        const double J3 = J3_i + J3_e;
        
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


        const int& meqn = source.getsize(2);

        if(_psi > meqn) continue;
        const double psi = qvals.get(i,_psi);
        source.set(i, _psi, -eps_Bcp_func(psi,time));

        if(_phi > meqn) continue;

        if( charge_density_coef )
        {
          const double phi = qvals.get(i,_phi);
          source.set(i, _phi, charge_density_coef*charge_density-eps_Ecp_func(phi,time));
        }
        else
        {
          source.set(i, _phi, 0.0);
        }

        if(_entropy_i > meqn) continue;

        source.set(i,_entropy_i, 0.0);
        source.set(i,_entropy_e, 0.0);

    }
*/

}
