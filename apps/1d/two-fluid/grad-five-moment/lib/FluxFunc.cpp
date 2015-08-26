#include "Components.h"         // Easier index into components of solution
#include "tensors.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     2D five-moment two-fluid equations
//
// The flux can be decoupled into a total of "three" smaller fluxes:
//
//      a) the flux for the electrons
//      b) the flux for the ions
//      c) the flux for Maxwell's equations
//
void FluxFunc(const dTensor1& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor2& flux)
{

    // First two blocks of flux function look like Euler (one for ions, and
    // one for electrons)
//  void FiveMomentFluxFunc( int n_offset, const dTensor2& Q, dTensor2& flux);
//  FiveMomentFluxFunc(0, Q, flux);     // Components for ions
//  FiveMomentFluxFunc(5, Q, flux);     // Components for electrons

//  // Last components describe Maxwell's equations
//  void MaxwellFluxFunc( int n_offset, const dTensor2& Q, dTensor2& flux);
//  MaxwellFluxFunc(10, Q, flux);

    // TODO - I'm not sure what this part is about ... -DS
//  void AdvectionFluxFunc( const dTensor2& Q, dTensor3& flux, int advIdx, int rhoIdx);
//  if(Q.getsize(2)<_entropy_i) return;
//  AdvectionFluxFunc(Q, flux, _entropy_i, _rho_i);
//  AdvectionFluxFunc(Q, flux, _entropy_e, _rho_e);

    // ---- This was copied directly from DoGPack ---- //
    // TODO - reuse FiveMomentFluxFunc ?? -DS

    const int numpts=xpts.getsize();
    double gamma,mass_ratio,debye,cs_light,larmor_radius;
    double rho_i,u1_i,u2_i,u3_i,press_i,energy_i;
    double rho_e,u1_e,u2_e,u3_e,press_e,energy_e;
    double B1,B2,B3,E1,E2,E3;

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);

        // Parameters
        const double mass_ratio     = global_ini_params.get_mass_ratio  ();
        const double debye          = global_ini_params.get_debye_length();
        const double larmor_radius  = global_ini_params.get_larmor_radius();
        const double gamma          = global_ini_params.get_gamma();
        const double cs_light       = global_ini_params.get_cs_light();

        // Variables
        rho_i    = Q.get(i,1);
        u1_i     = Q.get(i,2)/rho_i;
        u2_i     = Q.get(i,3)/rho_i;
        u3_i     = Q.get(i,4)/rho_i;
        energy_i = Q.get(i,5);

        rho_e    = Q.get(i,6);
        u1_e     = Q.get(i,7)/rho_e;
        u2_e     = Q.get(i,8)/rho_e;
        u3_e     = Q.get(i,9)/rho_e;
        energy_e = Q.get(i,10);

        B1       = Q.get(i,11);
        B2       = Q.get(i,12);
        B3       = Q.get(i,13);
        E1       = Q.get(i,14);
        E2       = Q.get(i,15);
        E3       = Q.get(i,16);

        press_i  = (gamma-1.0e0)*(energy_i - 0.5e0*rho_i*(u1_i*u1_i 
                    + u2_i*u2_i + u3_i*u3_i));

        press_e  = (gamma-1.0e0)*(energy_e - 0.5e0*rho_e*(u1_e*u1_e
                    + u2_e*u2_e + u3_e*u3_e));

        // Flux function
        flux.set(i,1,  rho_i*u1_i );
        flux.set(i,2,  rho_i*u1_i*u1_i + press_i );
        flux.set(i,3,  rho_i*u1_i*u2_i );
        flux.set(i,4,  rho_i*u1_i*u3_i );
        flux.set(i,5,  u1_i*(energy_i+press_i) ); 

        flux.set(i,6,  rho_e*u1_e );
        flux.set(i,7,  rho_e*u1_e*u1_e + press_e );
        flux.set(i,8,  rho_e*u1_e*u2_e );
        flux.set(i,9,  rho_e*u1_e*u3_e );
        flux.set(i,10, u1_e*(energy_e+press_e) );

        flux.set(i,11,  0.0 );
        flux.set(i,12, -E3 );
        flux.set(i,13,  E2 );
        flux.set(i,14,  0.0 );
        flux.set(i,15,  cs_light*cs_light*B3 );
        flux.set(i,16, -cs_light*cs_light*B2 );
    }


}

