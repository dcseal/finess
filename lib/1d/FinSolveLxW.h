#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

// ------------------------------------------------------------
// Function definitions
void SetBndValues(const dTensor2& node, dTensorBC2& aux, dTensorBC2& q);

void CopyQ(const dTensorBC2&, dTensorBC2&);

void ConSoln( 
    const dTensorBC2& aux,
    const dTensorBC2& q, 
    double t, string outputdir);

void BeforeStep(double dt, const dTensor2& node, dTensorBC2& aux, dTensorBC2& q);
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
//void SetRKinfo(int method2, RKinfo& rk);
//void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

// ------------------------------------------------------------
// LxW functions -- See LaxWendroff/
//  void ConstructIntegratedL( double dt, const dTensor2& node, 
//      const dTensorBC2& aux, const dTensorBC2& q,
//      dTensorBC2& f, dTensorBC2& fx, dTensorBC2& fxx, dTensorBC2&
//      qx, dTensorBC1& smax, dTensorBC2& F);
void ConstructIntegratedF( double dt, const dTensor2& node, 
    dTensorBC2& aux, dTensorBC2& q,
    dTensorBC1& smax, dTensorBC2& F);
void ConstructL(
        const dTensor2& node,
        const dTensorBC2& aux,
        const dTensorBC2& q,
        const dTensorBC2& F,  // <-- new term: integrated flux
        dTensorBC2& Lstar,
        dTensorBC1& smax);



//  void ConstructL(
//          const dTensor2& node,
//          dTensorBC2& aux,
//          dTensorBC2& q,      // setbndy conditions modifies q
//          dTensorBC2& Lstar,
//          dTensorBC1& smax);
// ------------------------------------------------------------

#endif
