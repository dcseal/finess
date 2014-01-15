#ifndef _LEGENDREPOLYS1D_H_
#define _LEGENDREPOLYS1D_H_
#include "dogdefs.h"
#include "dog_math.h"
#include "tensors.h"
#include "constants.h"

void evaluateLegendrePolys(const dTensor1& spts, dTensor2& phi );
void evaluateLegendrePolys( const double dx, const dTensor1& spts, dTensor2&
    phi, dTensor2& phi_x);
void setGaussLobattoPoints1d(dTensor1& x1d, dTensor1& wgt);
void setGaussLegendrePoints1d(dTensor1& x1d, dTensor1& w1d);


#endif
