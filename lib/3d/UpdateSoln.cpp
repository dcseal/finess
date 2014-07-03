#include "tensors.h"

// Update the solution using the constructed Lstar
//
// The Low-Storage RK methods use a combination of two different values of q,
// together with a right hand side, L.  This routine computes the following:
//
//     qnew = alpha1*qstar + alpha2*qnew + beta*dt*Lstar.
//
// The main loop covers all elements (including boundary cells) of the arrays.
//
// See also: RKinfo.h and SetRKinfo.cpp
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensorBC4& aux,
        const dTensorBC4& qstar, 
        const dTensorBC4& Lstar,
              dTensorBC4& qnew)
{

    const int numel = qnew.numel();
#pragma omp parallel for
    for( int k=0; k < numel; k++ )
    {
        double tmp = alpha1*qstar.vget(k)+alpha2*qnew.vget(k)+beta*dt*Lstar.vget(k);
        qnew.vset( k, tmp );
    }

    // Optional call to modify updated solution
//  void AfterUpdateSoln(const dTensor2& node,
//          const dTensorBC2& aux,
//          dTensorBC2& q,
//          double dt,
//          double beta);
//  AfterUpdateSoln(node, aux, qnew, dt, beta); 

}

// Update the solution using the constructed Lstar
//
// This version of UpdateSoln is used for the fifth-order time stepping.
void UpdateSoln(
    double g1,double g2, double g3, double delta, 
    double beta,double dt,
    const dTensorBC4& aux, const dTensorBC4& qold, const dTensorBC4& Lstar,
    dTensorBC4& q1, dTensorBC4& q2)
{

    const int     mx = q1.getsize(1);
    const int     my = q1.getsize(2);
    const int     mz = q1.getsize(3);
    const int   meqn = q1.getsize(4);
    const int   maux = aux.getsize(4);
    const int    mbc = q1.getmbc();

#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(my+mbc); j++)
    for (int k=(1-mbc); k<=(mz+mbc); k++)
    for (int m=1; m<=meqn; m++)        
    {

        double s1 =   q1.get(i,j,k,m);
        double s3 = qold.get(i,j,k,m);

        // update q2
        double s2 = q2.get(i,j,k,m) + delta*s1;
        q2.set(i,j,k,m, s2 );

        // update q
        double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.get(i,j,k,m);
        q1.set(i,j,k,m, tmp );
    }

    // Optional call to modify updated solution
//  void AfterUpdateSoln(const dTensor2& node,
//          const dTensorBC3& aux,
//          dTensorBC3& q,
//          double dt,
//          double beta);
//  AfterUpdateSoln(node,aux,q1,dt,beta); 
}
