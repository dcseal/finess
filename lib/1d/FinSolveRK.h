#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

// ------------------------------------------------------------
// Function definitions
void CopyQ(const dTensorBC2&, dTensorBC2&);

void ConSoln( const dTensor2& node, 
    const dTensorBC2& aux,
    const dTensorBC2& q, 
    double t, string outputdir);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, const dTensor2& node, dTensorBC2& aux, dTensorBC2& q);
void ConstructL( const dTensor2& node,
        dTensorBC2& aux,
        dTensorBC2& q,      // setbndy conditions modifies q
        dTensorBC2& Lstar,
        dTensorBC1& smax);
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensor2& node, 
        const dTensorBC2& aux,
        const dTensorBC2& qstar, 
        const dTensorBC2& Lstar,
        dTensorBC2& qnew);
void AfterStep(double dt, const dTensor2& node, dTensorBC2& aux, dTensorBC2& q);

double GetCFL(double dt, double dtmax,
        const dTensor1& prim_vol,
        const dTensorBC2& aux,
        const dTensorBC1& smax);

// TODO - add these calls in later
void BeforeFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);
void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);

// ------------------------------------------------------------
// Runge-Kutta information
void SetRKinfo(int method2, RKinfo& rk);
void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

#endif
