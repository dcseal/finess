#include <cmath>
#include "tensors.h"
#include "stdio.h"
#include "stdlib.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor1& xpts, 
	     dTensor2& auxvals)
{
  const int numpts=xpts.getsize();
  
  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);
      
      // bottom topography
      double bot = 0.5e0*exp(-100.0e0*pow((x-0.5),2));
      
      // derivative of bottom topography
      double dbdx = (100.0 - 200.0*x) * bot;
      
      //printf("x=%e,Auxvals=%e \n",x,dbdx);
      auxvals.set(i,1, dbdx );
    }
}
