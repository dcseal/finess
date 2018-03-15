#ifndef _CONSTRUCTS_H_
#define _CONSTRUCTS_H_

// This document is a compilation of all the "ConstructL" and variations of
// ConstructL.  The purpose of keeping it here is to consolidate all of the
// header files into a smaller document and make the code more legible.
//
// TODO - this doesn't actually have all of the construct files ... -DS 3/2018

#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "assert.h"

// Construct a time-integrated (averaged) right hand side function
void ConstructIntegratedR( double dt, const StateVars& Q,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1,
    const StateVars& Q1,
    double alpha2, double beta2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1, double charlie1, 
    const StateVars& Q1,
    double alpha2, double beta2, double charlie2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

// Construct a "local" time integration of the problem.  This is required to
// compute any time-average fluxes.
void LocalIntegrate( 
    int nterms, double dx, double dy, double xc, double yc,
    int meqn, int maux, int mpts_sten, int half_mpts_sten,
    const int i, const int j, const dTensorBC3& q, const dTensorBC3& aux, 
    const dTensorBC4& R, 
    dTensor1& f_t, dTensor1& f_tt,
    dTensor1& g_t, dTensor1& g_tt
    );

void ConstructLxWL( const StateVars& Q,
        const dTensorBC3& F, const dTensorBC3& G,
        dTensorBC3& Lstar, dTensorBC3& smax);

void ConstructLxWL( const StateVars& Q,
        const dTensorBC3& F, const dTensorBC3& G,
        dTensorBC3* pFhat, dTensorBC3* pGhat,
        dTensorBC3& Lstar, dTensorBC3& smax);

#endif
