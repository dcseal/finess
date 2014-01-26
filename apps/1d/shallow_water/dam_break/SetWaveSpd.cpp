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
    double hl,ul;
    double hr,ur;
    double h,u;
    
    // Left states
    hl = Ql.get(1);
    ul = Ql.get(2)/hl;
    
    // Right states
    hr = Qr.get(1);
    ur = Qr.get(2)/hr;
        
    // Average states
    h = 0.5e0*(hl + hr);
    u = 0.5e0*(ul + ur);
    
    // Minimum speed
    s1 = Min(u-sqrt(h),ul-sqrt(hl));
    
    // Maximum speed
    s2 = Max(u+sqrt(h),ur+sqrt(hr));
    
}
