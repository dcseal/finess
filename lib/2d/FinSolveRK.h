#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

// ------------------------------------------------------------
// Function definitions
void ConSoln( 
    const dTensorBC3& aux,
    const dTensorBC3& q, 
    double t, string outputdir);

double GetCFL(double dt, double dtmax,
        const dTensorBC3& aux,
        const dTensorBC2& smax);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, dTensorBC3& aux, dTensorBC3& q);
void ConstructL( 
        dTensorBC2& aux,
        dTensorBC2& q,      // setbndy conditions modifies q
        dTensorBC2& Lstar,
        dTensorBC1& smax);
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensorBC2& aux,
        const dTensorBC2& qstar, 
        const dTensorBC2& Lstar,
        dTensorBC2& qnew);
void AfterStep(double dt, dTensorBC2& aux, dTensorBC2& q);

// TODO - add these calls in later
void BeforeFullTimeStep(double dt, const dTensor1& prim_vol,
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);
void AfterFullTimeStep(double dt, const dTensor1& prim_vol,
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);

// ------------------------------------------------------------
// Runge-Kutta information
void SetRKinfo(int method2, RKinfo& rk);
void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

#endif
