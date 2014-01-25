#include "tensors.h"

// Update the solution using the constructed Lstar
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
        const dTensorBC3& aux,
        const dTensorBC3& qstar, 
        const dTensorBC3& Lstar,
              dTensorBC3& qnew)
{

    const int     mx = qnew.getsize(1);
    const int     my = qnew.getsize(2);
    const int   meqn = qnew.getsize(3);
    const int   maux =  aux.getsize(3);
    const int    mbc = qnew.getmbc();

#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(my+mbc); j++)
    for (int m=1; m<=meqn; m++)        
    {
        double tmp = alpha1*qstar.get(i,j,m) + alpha2*qnew.get(i,j,m)
            + beta*dt*Lstar.get(i,j,m);
        qnew.set(i, j, m, tmp );
    }

    // Optional call to modify updated solution
//  void AfterUpdateSoln(const dTensor2& node,
//          const dTensorBC2& aux,
//          dTensorBC2& q,
//          double dt,
//          double beta);
//  AfterUpdateSoln(node, aux, qnew, dt, beta); 

}
