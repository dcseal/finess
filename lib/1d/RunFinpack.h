#ifndef _RUN_FINESS_H_
#define _RUN_FINESS_H_

#include "StateVars.h"

void WriteQhelp( void );

void Output( const StateVars& Qstate, int nframe );

void QinitFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit( StateVars& Qstate );

void ConSoln( const StateVars& Qstate );

void FinSolveRK( StateVars& Qstate, double tend, double dtv[] );

void FinSolveLxW( StateVars& Qstate, double tend, double dtv[] );

void FinSolveMD( StateVars& Qstate, double tend, double dtv[] );

void FinSolveSDC( StateVars& Qstate, double tend, double dtv[] );

void FinSolveUser( StateVars& Qstate, double tend, double dtv[] );

// TODO - replace qin, auxin with StateVars ...
void SampleFunction( 
    int istart, int iend,
    const dTensorBC2& qin, 
    const dTensorBC2& auxin,  
    dTensorBC2& Fout,
    void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

// ------------------------------------------------------------

#endif
