#ifndef _RUN_FINESS_H_
#define _RUN_FINESS_H_

using namespace std;

// ------------------------------------------------------------
// Functions use in RunDogpack.cpp

void Output( const dTensorBC4& aux, const dTensorBC4& q,
        double t, int nframe, string outputdir);

void QinitFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit(dTensorBC4& aux, dTensorBC4& q);

void ConSoln(
    const dTensorBC4& aux, const dTensorBC4& q, 
    double t, string outputdir);

void FinSolveRK(
    dTensorBC4& aux, dTensorBC4& qold, dTensorBC4& qnew, 
    dTensorBC4& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir);

void FinSolveLxW(
    dTensorBC4& aux, dTensorBC4& qold, dTensorBC4& qnew, 
    dTensorBC4& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir);

void DogSolveUser(
        dTensorBC4& aux, dTensorBC4& qold, dTensorBC4& qnew,
        dTensorBC4& smax,
        double tstart, double tend,int nv, 
        double dtv[], const double cflv[],string outputdir);

void InitApp(IniDocument& ini_doc);

void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    int kstart, int kend,
    const dTensorBC4& qin, 
    const dTensorBC4& auxin,  
          dTensorBC4& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

// ------------------------------------------------------------

#endif
