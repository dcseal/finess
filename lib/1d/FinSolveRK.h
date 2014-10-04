#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

#include "StateVars.h"

// ------------------------------------------------------------
// Function definitions
void SetBndValues( StateVars& Q );
void CopyQ(const dTensorBC2&, dTensorBC2&);
void ConSoln( const StateVars& Qstate );

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, StateVars& Q );
void ConstructL( const StateVars& Q, dTensorBC2& Lstar, dTensorBC1& smax);
void ConstructL_NOC( const StateVars& Q, dTensorBC2& Lstar, dTensorBC1& smax);

// orders 1-4 time stepping:
void UpdateSoln(double alpha1, double alpha2, double beta, double dt, 
    const StateVars& Qstar, const dTensorBC2& Lstar, StateVars& Qnew );

// 5th-order time stepping (low-storage):
void UpdateSoln(
    double g1,double g2, double g3, double delta, 
    double beta,double dt,
    const StateVars& Qold,
    const dTensorBC2& Lstar,
    StateVars& Q1, StateVars& Q2);
void AfterStep(double dt, StateVars& Q );

double GetCFL(double dt, double dtmax,
        const dTensorBC2& aux,
        const dTensorBC1& smax);

void BeforeFullTimeStep(double dt, 
    const StateVars& Qold, StateVars& Qnew );
void AfterFullTimeStep(double dt, 
    const StateVars& Qold, StateVars& Qnew );

// ------------------------------------------------------------
// Runge-Kutta information
void SetRKinfo(int method2, RKinfo& rk);
void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

#endif
