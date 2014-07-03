#ifndef _CONSTRUCT_L_H_
#define _CONSTRUCT_L_H_

// --- User supplied functions --- //
void SetBndValues(dTensorBC4& aux, dTensorBC4& q);

void FluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor3& flux);
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, 
            const dTensor2& auxvals, dTensor2& source);

void ProjectLeftEig( int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
    const dTensor2& Qvals, dTensor2& Wvals);
void ProjectRightEig(int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
                     const dTensor2& Wvals, dTensor2& Qvals);

// Routines for computing maximum wave speeds

// Local wave speed
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge, 
    const dTensor1& Ql,   const dTensor1& Qr, 
    const dTensor1& Auxl, const dTensor1& Auxr,
    double& s1,double& s2);

// Global wave speed
void GlobalWaveSpd(
    const dTensorBC4& q, 
    const dTensorBC4& aux, 
    double& alpha1, double& alpha2, double& alpha3 );

// Routine for sampling a given function
void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    int kstart, int kend,
    const dTensorBC4& qin, 
    const dTensorBC4& auxin,  
          dTensorBC4& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

// Routine to deal with the silly mess where the Fluxes and the
// Projections are all defined separately.
void ConvertTranspose( const dTensor2& qin, dTensor2& qout );

#endif
