#include "dogdefs.h"
#include "StateVars.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH-ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void SetBndValues(StateVars& Q)
{

    dTensorBC4& q   = Q.ref_q();
    dTensorBC4& aux = Q.ref_aux();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int mz   = q.getsize(3);
    const int meqn = q.getsize(4);
    const int maux = aux.getsize(4);
    const int mbc  = q.getmbc();

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(0); i++)
    for (int j=1; j<=(my); j++)
    for (int k=1; k<=(mz); k++)
    {
	q.set(i, j, k, 1,  q.get(1-i, j, k, 1));
	q.set(i, j, k, 2, -q.get(1-i, j, k, 2));
	q.set(i, j, k, 3,  q.get(1-i, j, k, 3));
	q.set(i, j, k, 4,  q.get(1-i, j, k, 4));
	q.set(i, j, k, 5,  q.get(1-i, j, k, 5));
    }
    // ***********************************************


  // ***********************************************
  // RIGHT BOUNDARY
  // ***********************************************
#pragma omp parallel for
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=(my); j++)
    for (int k=1; k<=(mz); k++)    
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(mx,j,k,m);
            q.set(i,j,k,m, tmp );
        }
    }
    // ***********************************************


    // ***********************************************
    // FRONT BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(0); j++)    
    for (int k=1; k<=(mz); k++)
    {
	q.set(i, j, k, 1,  q.get(i, 1-j, k, 1));
	q.set(i, j, k, 2,  q.get(i, 1-j, k, 2));
	q.set(i, j, k, 3, -q.get(i, 1-j, k, 3));
	q.set(i, j, k, 4,  q.get(i, 1-j, k, 4));
	q.set(i, j, k, 5,  q.get(i, 1-j, k, 5));
    }
    // ***********************************************


    // ***********************************************
    // BACK BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(my); j<=(my+mbc); j++)    
    for (int k=1; k<=(mz); k++)
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(i,my,k,m);
            q.set(i,j,k,m, tmp );
        }	
        for (int m=1; m<=maux; m++)
        {
            double tmp = aux.get(i,my,k,m);
            aux.set(i,j,k,m, tmp );
        }	
    }
    // ***********************************************


    // ***********************************************
    // BOTTOM BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(my+mbc); j++)    
    for (int k=(1-mbc); k<=(0); k++)
    {
	q.set(i, j, k, 1,  q.get(i, j, 1-k, 1));
	q.set(i, j, k, 2,  q.get(i, j, 1-k, 2));
	q.set(i, j, k, 3,  q.get(i, j, 1-k, 3));
	q.set(i, j, k, 4, -q.get(i, j, 1-k, 4));
	q.set(i, j, k, 5,  q.get(i, j, 1-k, 5));
    }
    // ***********************************************
  

    // ***********************************************
    // TOP BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    for (int j=(1-mbc); j<=(my+mbc); j++)    
    for (int k=(mz+1); k<=(mz+mbc); k++)
    {
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(i,j,mz,m);
            q.set(i,j,k,m, tmp );
        }
        for (int m=1; m<=maux; m++)
        {
            double tmp = aux.get(i,j,mz,m);
            aux.set(i,j,k,m, tmp );
        }
    }
    // ***********************************************
}
