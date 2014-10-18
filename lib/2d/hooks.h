#ifndef _HOOKS_H_
#define _HOOKS_H_
#include "tensors.h"
#include "StateVars.h"

void AfterQinit( StateVars& Q );
void BeforeFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew );
void BeforeStep(double dt, StateVars& Q);
void AfterStep(double dt, StateVars& Q );
void AfterFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew );

#endif
