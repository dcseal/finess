#include "tensors.h"
#include "dog_math.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//     Set to zero so that f \pm \alpha q = f.
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
                const dTensor1& Auxr,
		double& s1,double& s2)
{
    // Minimum speed
//  s1 = Min(Auxl.get(1), 0.5e0*(Auxl.get(1)+Auxr.get(1)));
    
    // Maximum speed
//  s2 = Max(Auxr.get(1), 0.5e0*(Auxl.get(1)+Auxr.get(1)));

    s1 = 0.;
    s2 = 0.;

}
