#ifndef _CONSTRUCTL_H_
#define _CONSTRUCTL_H_

// User supplied functions:
void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
void ProjectLeftEig( const dTensor1& Aux_ave, const dTensor1& Q_ave, 
        const dTensor2& Qvals, dTensor2& Wvals);
void ProjectRightEig(const dTensor1& Aux_ave, const dTensor1& Q_ave, 
        const dTensor2& Wvals, dTensor2& Qvals);
void SetWaveSpd(const dTensor1& xedge, const dTensor1& Ql,
        const dTensor1& Qr, const dTensor1& Auxl, const dTensor1& Auxr,
        double& s1,double& s2);
void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
void SampleFunction( int istart, int iend,
        const dTensorBC2& qin, 
        const dTensorBC2& auxin,  dTensorBC2& Fout,
        void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

// Global wave speed
void GlobalWaveSpd( const dTensorBC2& q, const dTensorBC2& aux, double& alpha1 );

// Routine for querying the WENO reconstrution
void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);

// Routine to deal with the silly mess where the Fluxes and the
// Projections are all defined separately.
void ConvertTranspose( const dTensor2& qin, dTensor2& qout );

#endif
