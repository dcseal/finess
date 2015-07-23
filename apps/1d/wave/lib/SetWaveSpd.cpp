#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
// Two-fluid plasma  - TODO this is hacked together to look at speed of light
//                          only
//
void SetWaveSpd(const dTensor1& xedge,
        const dTensor1& Ql,
        const dTensor1& Qr,
        const dTensor1& Auxl,
        const dTensor1& Auxr,
        double& s1,double& s2)
{ 

    // We do not use this for the implicit solver
    s1 = global_ini_params.get_cs_light();
    s2 = s1;

}
