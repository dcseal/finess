#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

// ------------------------------------------------------------
// Function definitions
void ConSoln( const dTensorBC3& aux, const dTensorBC3& q, double t );

double GetCFL(double dt, double dtmax, const dTensorBC3& aux, const dTensorBC3& smax);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, dTensorBC3& aux, dTensorBC3& q);

// Used for orders 1--4:
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensorBC3& aux,
        const dTensorBC3& qstar, 
        const dTensorBC3& Lstar,
              dTensorBC3& qnew);

void SetBndValues(dTensorBC3& aux, dTensorBC3& q);

// Used for fifth-order stepper:
//  void UpdateSoln(
//      double g1,double g2, double g3, double delta, 
//      double beta,double dt,
//      const dTensorBC3& aux, const dTensorBC3& qold, const dTensorBC3& Lstar,
//      dTensorBC3& q1, dTensorBC3& q2);
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
// Lax-Wendroff information:
//
// TODO - copy over the ConstructIntegratedF over to this 
//        routine
//
// ------------------------------------------------------------

void ConstructLxWL(
        const dTensorBC3& aux,
        const dTensorBC3& q,
        dTensorBC3& F,         // <--- new term: integrated flux, f
        dTensorBC3& G,         // <--- new term: integrated flux, g
        dTensorBC3& Lstar,
        dTensorBC3& smax);

void ConstructIntegratedR( double dt, 
    const dTensorBC3& aux, const dTensorBC3& q,
    dTensorBC3& smax, 
    dTensorBC3& F, dTensorBC3& G);

#endif
