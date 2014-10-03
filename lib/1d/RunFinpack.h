#ifndef _RUN_FINESS_H_
#define _RUN_FINESS_H_

// ------------------------------------------------------------
// Functions use in RunDogpack.cpp
void Output(const dTensorBC2& aux, const dTensorBC2& q,
        double t, int nframe );

void QinitFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit(dTensorBC2& aux, dTensorBC2& q);

void ConSoln(
    const dTensorBC2& aux,
    const dTensorBC2& q, double t );

void FinSolveRK(
    dTensorBC2& aux, dTensorBC2& qnew, double tstart, 
    double tend, double dtv[] );

void FinSolveLxW(
    dTensorBC2& aux, dTensorBC2& qnew, double tstart, 
    double tend, double dtv[] );

void FinSolveMD(
    dTensorBC2& aux, dTensorBC2& qnew, double tstart, 
    double tend, double dtv[] );

void FinSolveSDC(
    dTensorBC2& aux, dTensorBC2& qnew, double tstart, 
    double tend, double dtv[] );

void FinSolveUser(
    dTensorBC2& aux, dTensorBC2& qnew, double tstart, 
    double tend, double dtv[] );

void SampleFunction( 
    int istart, int iend,
    const dTensorBC2& qin, 
    const dTensorBC2& auxin,  
    dTensorBC2& Fout,
    void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

// ------------------------------------------------------------

#endif
