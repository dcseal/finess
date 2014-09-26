#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

// ------------------------------------------------------------
// Function definitions
void SetBndValues(dTensorBC2& aux, dTensorBC2& q);

void CopyQ(const dTensorBC2&, dTensorBC2&);

void ConSoln( 
    const dTensorBC2& aux,
    const dTensorBC2& q, 
    double t, std::string outputdir);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, dTensorBC2& aux, dTensorBC2& q);
void ConstructL( const dTensorBC2& aux, const dTensorBC2& q,
        dTensorBC2& Lstar, dTensorBC1& smax);
void ConstructL_NOC( const dTensorBC2& aux, const dTensorBC2& q,
        dTensorBC2& Lstar, dTensorBC1& smax);

// orders 1-4 time stepping:
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensorBC2& aux,
        const dTensorBC2& qstar, 
        const dTensorBC2& Lstar,
        dTensorBC2& qnew);
// 5th-order time stepping (low-storage):
void UpdateSoln(
    double g1,double g2, double g3, double delta, 
    double beta,double dt,
    const dTensorBC2& aux,
    const dTensorBC2& qold, const dTensorBC2& Lstar,
    dTensorBC2& q1, dTensorBC2& q2);
void AfterStep(double dt, dTensorBC2& aux, dTensorBC2& q);

double GetCFL(double dt, double dtmax,
        const dTensorBC2& aux,
        const dTensorBC1& smax);

// TODO - add these calls in later
void BeforeFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);
void AfterFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q);

// ------------------------------------------------------------
// Runge-Kutta information
void SetRKinfo(int method2, RKinfo& rk);
void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

#endif
