#ifndef _RUN_FINESS_H_
#define _RUN_FINESS_H_

// ------------------------------------------------------------
// Functions use in RunDogpack.cpp
void Output(const dTensorBC2& aux, const dTensorBC2& q,
        double t, int nframe, std::string outputdir);

void QinitFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit(dTensorBC2& aux, dTensorBC2& q);

void ConSoln(
    const dTensorBC2& aux,
    const dTensorBC2& q, double t, std::string outputdir);

void FinSolveRK(
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], std::string outputdir);

void FinSolveLxW(
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], std::string outputdir);

void FinSolveMD(
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], std::string outputdir);

void FinSolveSDC(
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], std::string outputdir);

void DogSolveUser(
        dTensorBC2& aux, 
        dTensorBC2& qold,
        dTensorBC2& qnew,
        dTensorBC1& smax,
        double tstart, double tend,int nv, 
        double dtv[], const double cflv[], std::string outputdir);

void SampleFunction( 
    int istart, int iend,
    const dTensorBC2& qin, 
    const dTensorBC2& auxin,  
    dTensorBC2& Fout,
    void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

// ------------------------------------------------------------

#endif
