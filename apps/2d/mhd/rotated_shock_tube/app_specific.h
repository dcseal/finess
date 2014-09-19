#ifndef _APP_SPECIFIC_H_
#define _APP_SPECIFIC_H_
#include "tensors.h"
void AuxFunc(const dTensor2& xpts, 
	     dTensor2& auxvals);

void FluxFunc(const dTensor2& xpts,
	      const dTensor2& Q,
	      const dTensor2& Aux,
	      dTensor3& flux);


void QinitFunc(const dTensor2& xpts, dTensor2& qvals);
void SetBndValues( dTensorBC3& aux, dTensorBC3& q );
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals,
		    const dTensor2& auxvals, dTensor2& source);


#endif
