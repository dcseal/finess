#ifndef _FINSOLVE_SDC_H_
#define _FINSOLVE_SDC_H_

#include "StateVars.h"

void SetBndValues( StateVars& Q );
void CopyQ(const dTensorBC2&, dTensorBC2&);
void ConSoln( const StateVars& Qstate );

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

// SDC specific functions
void EulerStep(const double& dt, 
        dTensorBC1& smax, dTensorBC2& Lrhs,
        dTensorBC2& aux, dTensorBC2& qin, 
        dTensorBC2& qnew);

void TimeStepSDC(int method2, double t, double dt, dTensor1& dtvec, dTensor1& tvec);

void ResInt(double dt, 
        const dTensorBC2& L0, 
        const dTensorBC2& L1, 
        const dTensorBC2& L2, 
        const dTensorBC2& L3, 
        const dTensorBC2& L4, 
        const dTensorBC2& L5,
        dTensorBC3& ILout);

void StepSDCRK2(const double& dt, const int method[], const dTensor2& node,
        dTensorBC1& smax, dTensorBC3& Lrhs, dTensorBC3& Lstar,
        dTensorBC3& aux, dTensorBC3& qin, dTensorBC3& qstar,
        dTensorBC3& qnew);
void StepSDCdeltaRK2(const double& dt, const int method[], const dTensor2& node,
        dTensorBC1& smax, dTensorBC3& aux, dTensorBC3& qstar, dTensorBC3& Lstar,
        dTensorBC3& L1, dTensorBC3& L1new, dTensorBC3& L2, dTensorBC3& q1, 
        dTensorBC3& q2, int num, dTensorBC4& IL);

#endif
