#ifndef _APP_UTIL_H_
#define _APP_UTIL_H_

#include <cmath>

inline void ij_to_xy(double xlow, double ylow, double dx, double dy, int i, int j, double& x, double& y){
    x = xlow + (i - 0.5) * dx;
    y = ylow + (j - 0.5) * dy;
}


void UpdateSolnWithAux(double alpha1, double alpha2, double beta, double dt,
    const StateVars& Qstar, const dTensorBC3& Lstar, const dTensorBC3& Lauxstar, StateVars& Qnew);


#endif
