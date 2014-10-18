#ifndef _APP_SPECIFIC_H_
#define _APP_SPECIFIC_H_

#include "tensors.h"
#include "StateVars.h"

void QinitFunc(const dTensor2& xpts, dTensor2& qvals);
void AuxFunc(const dTensor2& xpts, dTensor2& auxvals);
void FluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor3& flux);
void SetBndValues( StateVars& Q );
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, const dTensor2& auxvals, dTensor2& source);

#endif
