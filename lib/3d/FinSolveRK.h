#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

// ------------------------------------------------------------
// Function definitions
void ConSoln( 
    const dTensorBC4& aux,
    const dTensorBC4& q, 
    double t, string outputdir);

double GetCFL(double dt, double dtmax,
        const dTensorBC4& aux,
        const dTensorBC4& smax);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q);
void ConstructL( 
        dTensorBC4& aux,
        dTensorBC4& q,      // setbndy conditions modifies q
        dTensorBC4& Lstar,
        dTensorBC4& smax);

// Used for orders 1--4:
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensorBC4& aux,
        const dTensorBC4& qstar, 
        const dTensorBC4& Lstar,
              dTensorBC4& qnew);

// Used for fifth-order stepper:
void UpdateSoln(
    double g1,double g2, double g3, double delta, 
    double beta,double dt,
    const dTensorBC4& aux, const dTensorBC4& qold, const dTensorBC4& Lstar,
    dTensorBC4& q1, dTensorBC4& q2);
void AfterStep(double dt, dTensorBC4& aux, dTensorBC4& q);

// Called once before each full time step (c.f. BeforeStep, which is called
// before each stage in an RK method)
void BeforeFullTimeStep(double dt, 
               dTensorBC4& auxold, dTensorBC4& aux, 
               dTensorBC4& qold,   dTensorBC4& q);

void AfterFullTimeStep(double dt, 
               dTensorBC4& auxold, dTensorBC4& aux, 
               dTensorBC4& qold,   dTensorBC4& q);

// ------------------------------------------------------------
// Runge-Kutta information
void SetRKinfo(int method2, RKinfo& rk);
void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

#endif
