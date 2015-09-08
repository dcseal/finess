#ifndef _RUN_FINESS_H_
#define _RUN_FINESS_H_

#include "StateVars.h"

void WriteQhelp( void );

void Output( const StateVars& Qnew, int nframe );
void QinitFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit( StateVars& Qnew );

void ConSoln( const StateVars& Q );

// Time stepping methods:
void FinSolveRK     ( StateVars& Qnew, double tend, double dtv[] );
void FinSolveLxW    ( StateVars& Qnew, double tend, double dtv[] );
void FinSolveMD     ( StateVars& Qnew, double tend, double dtv[] );
void FinSolveUser   ( StateVars& Qnew, double tend, double dtv[] );

void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC3& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

// ------------------------------------------------------------

#endif
