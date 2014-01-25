#ifndef _RUN_FINESS_H_
#define _RUN_FINESS_H_

using namespace std;

// ------------------------------------------------------------
// Functions use in RunDogpack.cpp
void GridSetup( dTensor2& node, dTensor1& prim_vol);

void Output(const dTensor2& node, const dTensorBC2& aux, const dTensorBC2& q,
        double t, int nframe, string outputdir);

void QinitFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit(const dTensor2& node, dTensorBC2& aux, dTensorBC2& q);

void ConSoln(const dTensor2& node, 
    const dTensorBC2& aux,
    const dTensorBC2& q, double t, string outputdir);

void FinSolveRK(
    const dTensor2& node, const dTensor1& prim_vol, 
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir);

void DogSolveUser(const dTensor2& node, 
        const dTensor1& prim_vol, 
        dTensorBC2& aux, 
        dTensorBC2& qold,
        dTensorBC2& qnew,
        dTensorBC1& smax,
        double tstart, double tend,int nv, 
        double dtv[], const double cflv[],string outputdir);

void InitApp(IniDocument& ini_doc);

void SampleFunction( 
    int istart, int iend,
    const dTensor2& node,
    const dTensorBC2& qin, 
    const dTensorBC2& auxin,  
    dTensorBC2& Fout,
    void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

// ------------------------------------------------------------

#endif
