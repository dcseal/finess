#ifndef _APP_UTIL_H_
#define _APP_UTIL_H_

#include <cmath>


inline void ij_to_xy(double xlow, double ylow, double dx, double dy, int i, int j, double& x, double& y){
    x = xlow + (i - 0.5) * dx;
    y = ylow + (j - 0.5) * dy;
}

inline double A3_exact(double angle, double t, double x, double y){
    using std::sin;
    using std::cos;
    return -x * sin(angle) + y * cos(angle) + 0.1 / (2*M_PI) * cos(2*M_PI *(x*cos(angle) + y*sin(angle) + t));
}

void UpdateSolnWithAux(double alpha1, double alpha2, double beta, double dt,
    const StateVars& Qstar, const dTensorBC3& Lstar, const dTensorBC3& Lauxstar, StateVars& Qnew);


#endif
