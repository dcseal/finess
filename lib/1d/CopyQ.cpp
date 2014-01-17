// TODO - There are options in the Tensor class that can make this entire call
// obsolete ...

#include <iostream>
#include "tensors.h"
using namespace std;

// Copy contents of qin into qout
void CopyQ(const dTensorBC2& qin, dTensorBC2& qout)
{
    const int mbc = qin.getmbc();
    const int imax = qin.getsize(1);
    const int mmax = qin.getsize(2);

#pragma omp parallel for
    for (int i=1-mbc; i<=(imax+mbc); i++)    
    for (int m=1; m<=mmax; m++)        
	{
	  double tmp1 = qin.get(i,m);
	  qout.set(i,m, tmp1);
	}
    
}

// Copy contents of qin into qout
void CopyQ(const dTensorBC3& qin, dTensorBC3& qout)
{
  const int mbc = qin.getmbc();
  const int imax = qin.getsize(1);
  const int mmax = qin.getsize(2);
  const int kmax = qin.getsize(3);
  
#pragma omp parallel for
  for (int i=1-mbc; i<=(imax+mbc); i++)    
    for (int m=1; m<=mmax; m++)        
      for (int k=1; k<=kmax; k++)
	{
	  double tmp1 = qin.get(i,m,k);
	  qout.set(i,m,k, tmp1);
	}
    
}

// Copy contents of qin(:,:,:,mopt) into qout(:,:,:).  
void CopyQ(const int& mopt, const dTensorBC4& qin, dTensorBC3& qout)
{
  const int mbc = qin.getmbc();
  const int imax = qin.getsize(1);
  const int mmax = qin.getsize(2);
  const int kmax = qin.getsize(3);

#pragma omp parallel for
  for (int i=1-mbc; i<=(imax+mbc); i++)
    for (int m=1; m<=mmax; m++)
      for (int k=1; k<=kmax; k++)
	{
	  double tmp1 = qin.get(i,m,k,mopt);
	  qout.set(i,m,k, tmp1);
	}    
}

// Copy contents of qin(:,:,:) into qout(:,:,:,mopt).  
void CopyQ(const int& mopt, const dTensorBC3& qin, dTensorBC4& qout)
{    
  const int mbc = qin.getmbc();
  const int imax = qin.getsize(1);
  const int mmax = qin.getsize(2);
  const int kmax = qin.getsize(3);

#pragma omp parallel for
  for (int i=1-mbc; i<=(imax+mbc); i++)
    for (int m=1; m<=mmax; m++)
      for (int k=1; k<=kmax; k++)
	{
	  double tmp1 = qin.get(i,m,k);
	  qout.set(i,m,k, mopt, tmp1);
	}

}
