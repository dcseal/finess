#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

#include <string>

#include "StateVars.h"

// ------------------------------------------------------------
// Function definitions
void ConSoln( const StateVars& Q );

double GetCFL(double dt, double dtmax,
        const dTensorBC3& aux,
        const dTensorBC3& smax);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, StateVars& Q );
void ConstructL( StateVars& Q, dTensorBC3& Lstar, dTensorBC3& smax);

// Used for orders 1--4:
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
    const StateVars& Q1, const dTensorBC3& Lstar, StateVars& Qnew);

// Used for fifth-order stepper:
void UpdateSoln(
    double g1, double g2, double g3, double delta, 
    double beta, double dt,
    const StateVars& Qold, const dTensorBC3& Lstar,
    StateVars& Q1, StateVars& Q2);

void AfterStep(double dt, StateVars& Q );

// Called once before each full time step (c.f. BeforeStep, which is called
// before each stage in an RK method)
void BeforeFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew);
void AfterFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew);

// ------------------------------------------------------------
// Runge-Kutta information
void SetRKinfo(int method2, RKinfo& rk);
void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

#endif
