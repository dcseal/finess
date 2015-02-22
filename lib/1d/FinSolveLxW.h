#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

#include "StateVars.h"

void SetBndValues( StateVars& Q );
void ConSoln( const StateVars& Qstate );

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, StateVars& Q );
void ConstructL( const StateVars& Q, dTensorBC2& Lstar, dTensorBC1& smax);
void AfterStep(double dt, StateVars& Q );

double GetCFL(double dt, double dtmax,
        const dTensorBC2& aux,
        const dTensorBC1& smax);

void BeforeFullTimeStep(double dt, 
    const StateVars& Qold, StateVars& Qnew );
void AfterFullTimeStep(double dt, 
    const StateVars& Qold, StateVars& Qnew );


// ------------------------------------------------------------
// LxW functions -- See LaxWendroff/
void ConstructIntegratedF( double dt, const StateVars& Q, dTensorBC1& smax, dTensorBC2& F);
void ConstructLxWL(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        const dTensorBC2& F,  // <-- new term: integrated flux
        dTensorBC2& Lstar,
        dTensorBC1& smax);

void ConstructLxWL(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        const dTensorBC2& F,    
        dTensorBC2& fhat,
        dTensorBC2& Lstar,
        dTensorBC1& smax);

void ConstructLFL( const StateVars& Q, dTensorBC2& fhat );

void ApplyMPPLimiter1D( const double dt, const dTensorBC2& q,
    const dTensorBC2& fLF, dTensorBC2& fHat );

#endif
