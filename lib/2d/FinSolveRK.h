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
        dTensorBC3& aux,
        dTensorBC3& q,      // setbndy conditions modifies q
        dTensorBC3& Lstar,
        dTensorBC2& smax);
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensorBC3& aux,
        const dTensorBC3& qstar, 
        const dTensorBC3& Lstar,
              dTensorBC3& qnew);
void AfterStep(double dt, dTensorBC3& aux, dTensorBC3& q);

// Called once before each full time step (c.f. BeforeStep, which is called
// before each stage in an RK method)
void BeforeFullTimeStep(double dt, 
               dTensorBC3& auxold, dTensorBC3& aux, 
               dTensorBC3& qold,   dTensorBC3& q);

void AfterFullTimeStep(double dt, 
               dTensorBC3& auxold, dTensorBC3& aux, 
               dTensorBC3& qold,   dTensorBC3& q);

// ------------------------------------------------------------
// Runge-Kutta information
void SetRKinfo(int method2, RKinfo& rk);
void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

#endif
