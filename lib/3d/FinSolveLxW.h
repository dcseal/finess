#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

#include "StateVars.h"

void ConSoln( const StateVars& Q );
double GetCFL(double dt, double dtmax,
        const dTensorBC4& aux,
        const dTensorBC4& smax);

// ------------------------------------------------------------
// Lax-Wendroff information:
//
// TODO - copy over the ConstructIntegratedF over to this 
//        routine
//
// ------------------------------------------------------------

void ConstructIntegratedR( double dt, const StateVars& Q, 
    dTensorBC4& smax, dTensorBC4& F, dTensorBC4& G, dTensorBC4& H);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, StateVars& Q );
void SetBndValues(StateVars& Q );
void AfterStep(double dt, StateVars& Q);

// Called once before each full time step (c.f. BeforeStep, which is called
// before each stage in an RK method)
void BeforeFullTimeStep(double dt, StateVars& Qold, StateVars& Q );
void AfterFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew );
void ConstructLFL( const double dt, const StateVars& Qnew, 
    dTensorBC4& fLF, dTensorBC4& gLF, dTensorBC4& hLF, dTensorBC4& Lstar, dTensorBC4& smax );
void ApplyMPPLimiter3D( 
        const double dt, const dTensorBC4& q, 
        const dTensorBC4& fLF, const dTensorBC4& gLF, const dTensorBC4& hLF,
        dTensorBC4& fHat, dTensorBC4& gHat, dTensorBC4& hHat);

void ConstructLxWL( const StateVars& Q,
        const dTensorBC4& F, const dTensorBC4& G, const dTensorBC4& H,
        dTensorBC4& Lstar, dTensorBC4& smax);

void ConstructLxWL( const StateVars& Q,
        const dTensorBC4& F,         // <--- new term: integrated flux, f
        const dTensorBC4& G,         // <--- new term: integrated flux, g
        const dTensorBC4& H,
	dTensorBC4* pFhat, dTensorBC4* pGhat, dTensorBC4* pHhat,
        dTensorBC4& Lstar,
        dTensorBC4& smax);

#endif
