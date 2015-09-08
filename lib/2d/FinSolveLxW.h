#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

#include "StateVars.h"

// ------------------------------------------------------------
// Function definitions
void ConSoln( const StateVars& Q );

double GetCFL(double dt, double dtmax, const dTensorBC3& aux, const dTensorBC3& smax);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, StateVars& Q );
void SetBndValues(StateVars& Q );
void AfterStep(double dt, StateVars& Q);

// Called once before each full time step (c.f. BeforeStep, which is called
// before each stage in an RK method)
void BeforeFullTimeStep(double dt, StateVars& Qold, StateVars& Q );
void AfterFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew );

// ------------------------------------------------------------
// Lax-Wendroff information:
//
// TODO - copy over the ConstructIntegratedF over to this 
//        routine
//
// ------------------------------------------------------------

void ConstructLxWL( const StateVars& Q,
        const dTensorBC3& F,         // <--- new term: integrated flux, f
        const dTensorBC3& G,         // <--- new term: integrated flux, g
        dTensorBC3& Lstar,
        dTensorBC3& smax);

void ConstructLxWL( const StateVars& Q,
        dTensorBC3& F,     dTensorBC3& G,
        dTensorBC3* pFhat, dTensorBC3* pGhat,
        dTensorBC3& Lstar, dTensorBC3& smax);

void ConstructIntegratedR( double dt, const StateVars& Q, 
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

// Construct a Lax-Friedrich's flux
void ConstructLFL( const double dt, StateVars& Q,
    dTensorBC3& Fhat, dTensorBC3& Ghat,
    dTensorBC3& Lstar, dTensorBC3& smax);

void ApplyMPPLimiter2D( 
        const double dt, const dTensorBC3& q, 
        const dTensorBC3& fLF, const dTensorBC3& gLF,
        dTensorBC3& fHat, dTensorBC3& gHat );

#endif
