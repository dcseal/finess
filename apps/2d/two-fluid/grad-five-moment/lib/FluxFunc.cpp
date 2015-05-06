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
void FluxFunc(const dTensor2& xpts,
	      const dTensor2& Q,
	      const dTensor2& Aux,
	      dTensor3& flux)
{

    // First two blocks of flux function look like Euler (one for ions, and
    // one for electrons)
    void FiveMomentFluxFunc( int n_offset, const dTensor2& Q, dTensor3& flux);
    FiveMomentFluxFunc(0, Q, flux);     // Components for ions
    FiveMomentFluxFunc(5, Q, flux);     // Components for electrons

    // Last components describe Maxwell's equations
    void MaxwellFluxFunc( int n_offset, const dTensor2& Q, dTensor3& flux);
    MaxwellFluxFunc(10, Q, flux);

    // TODO - I'm not sure what this part is about ... -DS
//  void AdvectionFluxFunc( const dTensor2& Q, dTensor3& flux, int advIdx, int rhoIdx);
//  if(Q.getsize(2)<_entropy_i) return;
//  AdvectionFluxFunc(Q, flux, _entropy_i, _rho_i);
//  AdvectionFluxFunc(Q, flux, _entropy_e, _rho_e);

}

void FluxFunc1( const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor2& flux)
{

    void FiveMomentFluxFunc1( int n_offset, const dTensor2& Q, dTensor2& flux);
    FiveMomentFluxFunc1(0, Q, flux);
    FiveMomentFluxFunc1(5, Q, flux);

    void MaxwellFluxFunc1( int n_offset, const dTensor2& Q, dTensor2& flux);
    MaxwellFluxFunc1(10, Q, flux);

    // TODO - I'm not sure what this part is about ... -DS
//  void AdvectionFluxFunc1( const dTensor2& Q, dTensor2& flux, int advIdx, int rhoIdx);
//  if(Q.getsize(2)<_entropy_i) return;
//  AdvectionFluxFunc1(Q, flux, _entropy_i, _rho_i);
//  if(Q.getsize(2)<_entropy_e) return;
//  AdvectionFluxFunc1(Q, flux, _entropy_e, _rho_e);

}

void FluxFunc2( const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor2& flux)
{

    void FiveMomentFluxFunc2( int n_offset, const dTensor2& Q, dTensor2& flux);
    FiveMomentFluxFunc2(0, Q, flux);
    FiveMomentFluxFunc2(5, Q, flux);

    void MaxwellFluxFunc2( int n_offset, const dTensor2& Q, dTensor2& flux);
    MaxwellFluxFunc2(10, Q, flux);

    // TODO - I'm not sure what this part is about ... -DS
//  void AdvectionFluxFunc2( const dTensor2& Q, dTensor2& flux, int advIdx, int rhoIdx);
//  if(Q.getsize(2)<_entropy_i) return;
//  AdvectionFluxFunc2(Q, flux, _entropy_i, _rho_i);
//  if(Q.getsize(2)<_entropy_e) return;
//  AdvectionFluxFunc2(Q, flux, _entropy_e, _rho_e);

}

