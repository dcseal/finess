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

void ConstructLxWL( const StateVars& Q,
        const dTensorBC4& F, const dTensorBC4& G, const dTensorBC4& H,
        dTensorBC4& Lstar, dTensorBC4& smax);

#endif
