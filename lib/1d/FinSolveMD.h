#ifndef __FINSOLVE_MD__
#define __FINSOLVE_MD__

#include "StateVars.h"

void SetBndValues( StateVars& Q );
void ConSoln( const StateVars& Qstate );

void BeforeStep(double dt, StateVars& Q );
void ConstructLxWL(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        const dTensorBC2& F,  // <-- new term: integrated flux
        dTensorBC2& Lstar,
        dTensorBC1& smax);
void AfterStep(double dt, StateVars& Q );

double GetCFL(double dt, double dtmax,
        const dTensorBC2& aux,
        const dTensorBC1& smax);

void BeforeFullTimeStep(double dt, 
    const StateVars& Qold, StateVars& Qnew );
void AfterFullTimeStep(double dt, 
    const StateVars& Qold, StateVars& Qnew );


// ------------------------------------------------------------
// Taylor series integration
void ConstructIntegratedF( double dt, 
    dTensorBC2& aux, dTensorBC2& q,
    dTensorBC1& smax, dTensorBC2& F);
// ------------------------------------------------------------

// ------------------------------------------------------------
// Multiderivative integration
//
// These functions are for the two-stage methods.  One contains
// two-derivatives, and the second contains three derivatives.
void ConstructIntegratedF( double dt, 
    double alpha1, double beta1,
    dTensorBC2& aux1, dTensorBC2& q1,
    double alpha2, double beta2,
    dTensorBC2& aux2, dTensorBC2& q2,
    dTensorBC1& smax, dTensorBC2& F);

void ConstructIntegratedF( double dt, 
    double alpha1, double beta1, double charlie1,
    dTensorBC2& aux1, dTensorBC2& q1,
    double alpha2, double beta2, double charlie2,
    dTensorBC2& aux2, dTensorBC2& q2,
    dTensorBC1& smax, dTensorBC2& F);
// ------------------------------------------------------------

#endif 
