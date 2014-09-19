#ifndef _HOOKS_H_
#define _HOOKS_H_
#include "tensors.h"
void AfterQinit( dTensorBC3& aux, dTensorBC3& q );
void BeforeFullTimeStep(double dt, 
		       dTensorBC3& auxold, dTensorBC3& aux, 
		       dTensorBC3& qold,   dTensorBC3& q);

void BeforeStep(double dt, dTensorBC3& aux, dTensorBC3& q);
void AfterStep(double dt, dTensorBC3& aux, dTensorBC3& q);
void AfterFullTimeStep(double dt,
		       dTensorBC3& auxold, dTensorBC3& aux, 
		       dTensorBC3& qold,   dTensorBC3& q);
#endif
