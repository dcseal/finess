#ifndef _CONSTRUCTHJ_L_H_
#define _CONSTRUCTHJ_L_H_

#include "StateVars.h"

void ConstructHJ_L_Order1(const StateVars& Q, dTensorBC3& Lauxstar);
void ConstructHJ_L_Order3(const StateVars& Q, dTensorBC3& Lauxstar, double dt);


#endif
