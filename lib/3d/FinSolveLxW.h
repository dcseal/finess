#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

double GetCFL(double dt, double dtmax,
        const dTensorBC4& aux,
        const dTensorBC4& smax);

// ------------------------------------------------------------
// Lax-Wendroff information:
//
// TODO - copy over the ConstructIntegratedF over to this 
//        routine
//
// ------------------------------------------------------------

void ConstructIntegratedR( double dt, 
    const dTensorBC3& aux, const dTensorBC3& q,
    dTensorBC3& smax, 
    dTensorBC3& F, dTensorBC3& G);

#endif
