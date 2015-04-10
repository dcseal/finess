#include "dog_math.h"
#include "tensors.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Simple advection equation
//
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        double& s1,double& s2)
{

    // Normal component of velocity
    double un_left  = nvec.get(1)*Ql.get(1) + nvec.get(2)*Ql.get(1) + nvec.get(3)*Ql.get(1);
    double un_right = nvec.get(1)*Qr.get(1) + nvec.get(2)*Qr.get(1) + nvec.get(3)*Qr.get(1);
    double un_av    = 0.5*(un_left + un_right);

    // Minimum speed
    s1 = Min(un_left, un_av);

    // Maximum speed
    s2 = Max(un_right, un_av);

}
