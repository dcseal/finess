#include <cmath>
#include "dog_math.h"
#include "tensors.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Shallow water equations
//
void SetWaveSpd(const dTensor1& xedge,
        const dTensor1& Ql,
        const dTensor1& Qr,
        const dTensor1& Auxl,
        const dTensor1& Auxr,
        double& s1,double& s2)
{

    // Left states
    double hl = Ql.get(1);
    double ul = Ql.get(2)/hl;

    // Right states
    double hr = Qr.get(1);
    double ur = Qr.get(2)/hr;

    // Average states
    double h = 0.5e0*(hl + hr);
    double u = 0.5e0*(ul + ur);

    // Minimum speed
    s1 = Min(u-sqrt(h),ul-sqrt(hl));

    // Maximum speed
    s2 = Max(u+sqrt(h),ur+sqrt(hr));

}
