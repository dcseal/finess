#include "tensors.h"

// Update the solution using the constructed Lstar
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensorBC2& aux,
        const dTensorBC2& qstar, 
        const dTensorBC2& Lstar,
        dTensorBC2& qnew)
{
    const int melems = qnew.getsize(1);
    const int   meqn = qnew.getsize(2);
    const int   maux = aux.getsize(2);
    const int mbc = qnew.getmbc();

#pragma omp parallel for
    for (int j=(1-mbc); j<=(melems+mbc); j++)
    for (int m=1; m<=meqn; m++)        
    {
        double tmp = alpha1*qstar.get(j,m) + alpha2*qnew.get(j,m)
            + beta*dt*Lstar.get(j,m);
        qnew.set(j, m, tmp );
    }

    // Optional call to modify updated solution
//  void AfterUpdateSoln(
//          const dTensorBC2& aux,
//          dTensorBC2& q,
//          double dt,
//          double beta);
//  AfterUpdateSoln( aux, qnew, dt, beta); 

}

// Update the solution using the constructed Lstar
void UpdateSoln(
    double g1,double g2, double g3, double delta, 
    double beta,double dt,
    const dTensorBC2& aux,
    const dTensorBC2& qold, const dTensorBC2& Lstar,
    dTensorBC2& q1, dTensorBC2& q2)
{
    const int numel = q1.numel();
#pragma omp parallel for
    for( int k=0; k < numel; k++ )
    {
        double s1 = q1.vget( k );
        double s3 = qold.vget( k );

        // update q2
        double s2 = q2.vget( k ) + delta*s1;
        q2.vset(k, s2 );

        // update q
        double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.vget(k);
        q1.vset(k, tmp );

    }

    // Optional call to modify updated solution
    //
    // TODO - include this call
    //
//  void AfterUpdateSoln(
//          const dTensorBC3& aux,
//          dTensorBC3& q,
//          double dt,
//          double beta);
//  AfterUpdateSoln(aux,q1,dt,beta); 

}
