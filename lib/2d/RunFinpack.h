#ifndef _RUN_FINESS_H_
#define _RUN_FINESS_H_

using namespace std;

// ------------------------------------------------------------
// Functions use in RunDogpack.cpp

void Output( const dTensorBC3& aux, const dTensorBC3& q,
        double t, int nframe, string outputdir);

void QinitFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit(dTensorBC3& aux, dTensorBC3& q);

void ConSoln(
    const dTensorBC3& aux, const dTensorBC3& q, 
    double t, string outputdir);

void FinSolveRK(
    dTensorBC3& aux, dTensorBC3& qold, dTensorBC3& qnew, 
    dTensorBC3& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir);

void FinSolveLxW(
    dTensorBC3& aux, dTensorBC3& qold, dTensorBC3& qnew, 
    dTensorBC3& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir);


void DogSolveUser(
        dTensorBC3& aux, dTensorBC3& qold, dTensorBC3& qnew,
        dTensorBC3& smax,
        double tstart, double tend,int nv, 
        double dtv[], const double cflv[],string outputdir);

void InitApp(IniDocument& ini_doc);

void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC3& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

// ------------------------------------------------------------

#endif
