#ifndef _MISC2D_H_
#define _MISC2D_H_

#include "tensors.h"
#include "StateVars.h"

double GetCFL(double dt, double dtmax, const dTensorBC3& aux, const dTensorBC3& smax);
double GetCFL(double dt, double dtmax, double alpha1, double alpha2 );
void ConSoln( const StateVars& Q );

#endif
