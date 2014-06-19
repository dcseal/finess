#ifndef _GLOBAL_WAVE_SPD_H_
#define _GLOBAL_WAVE_SPD_H_

// Compute a global wave speed.

// Note: each Numerical flux should be consistent.  That is, if Ql = Qr, then
// fhat( Ql, Qr ) = f( Ql ), where f is the real flux.  For now, we will rely on
// this fact in order to construct the correct maximum wave speeds.
void SetWaveSpd(const dTensor1& nvec, 
        const dTensor1& xedge,
        const dTensor1& Ql, 
        const dTensor1& Qr,
        const dTensor1& Auxl, 
        const dTensor1& Auxr,
        double& s1,double& s2);

#endif
