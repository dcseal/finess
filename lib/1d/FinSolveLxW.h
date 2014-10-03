#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

#include "StateVars.h"

// ------------------------------------------------------------
// Function definitions
void SetBndValues(dTensorBC2& aux, dTensorBC2& q);

void CopyQ(const dTensorBC2&, dTensorBC2&);

void ConSoln( const StateVars& Qstate );

void BeforeStep(double dt, dTensorBC2& aux, dTensorBC2& q);
void AfterStep(double dt, dTensorBC2& aux, dTensorBC2& q);

double GetCFL(double dt, double dtmax,
        const dTensorBC2& aux,
        const dTensorBC1& smax);

// TODO - add these calls in later
void BeforeFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);
void AfterFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);

// ------------------------------------------------------------
// Runge-Kutta information
//void SetRKinfo(int method2, RKinfo& rk);
//void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

// ------------------------------------------------------------
// LxW functions -- See LaxWendroff/
void ConstructIntegratedF( double dt, 
    dTensorBC2& aux, dTensorBC2& q,
    dTensorBC1& smax, dTensorBC2& F);
void ConstructLxWL(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        const dTensorBC2& F,  // <-- new term: integrated flux
        dTensorBC2& Lstar,
        dTensorBC1& smax);


#endif
