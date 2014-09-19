#ifndef _CONSTRUCTS_H_
#define _CONSTRUCTS_H_

#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "assert.h"

void ConstructIntegratedR( double dt, 
    const dTensorBC3& aux, const dTensorBC3& q,
    dTensorBC3& smax, 
    dTensorBC3& F, dTensorBC3& G);


void LocalIntegrate( 
    int nterms, double dx, double dy, double xc, double yc,
    int meqn, int maux, int mpts_sten, int half_mpts_sten,
    const int i, const int j, const dTensorBC3& q, const dTensorBC3& aux, 
    const dTensorBC4& R, 
    dTensor1& f_t, dTensor1& f_tt,
    dTensor1& g_t, dTensor1& g_tt
    );

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1,
    const dTensorBC3& aux1, const dTensorBC3& q1,
    double alpha2, double beta2,
    const dTensorBC3& aux2, const dTensorBC3& q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);


void ConstructIntegratedR( double dt, 
    double alpha1, double beta1, double charlie1, 
    const dTensorBC3& aux1, const dTensorBC3& q1,
    double alpha2, double beta2, double charlie2,
    const dTensorBC3& aux2, const dTensorBC3& q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);
void ConstructLxWL(
        const dTensorBC3& aux,
        const dTensorBC3& q,
        dTensorBC3& F,         // <--- new term: integrated flux, f
        dTensorBC3& G,         // <--- new term: integrated flux, g
        dTensorBC3& Lstar,
        dTensorBC3& smax);
#endif
