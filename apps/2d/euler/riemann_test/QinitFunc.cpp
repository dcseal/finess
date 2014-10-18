#include <fstream>
#include "dogdefs.h"
#include "IniParams.h"
#include "dog_math.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
  const int numpts=xpts.getsize(1);
  const double gamma = global_ini_params.get_gamma();
  const int OPT = global_ini_params.get_OPT();

  // Loop over grid points
  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i,1);
      double y = xpts.get(i,2);
	
      double rho,u1,u2,u3,press;
      if (OPT==1)
	{
	  if(x>0.5)
	    {
	      rho   =  1.0;
	      u1    =  0.0;
	      u2    =  0.0;
	      u3    =  0.0;
	      press =  1.0;
	    }
	  else
	    {
	      rho   =  3.0;
	      u1    =  0.0;
	      u2    =  0.0;
	      u3    =  0.0;
	      press =  3.0;
	    } 
	}
      else
	{
	  if(y>0.5)
            {
	      rho   =  1.0;
	      u1    =  0.0;
	      u2    =  0.0;
	      u3    =  0.0;
	      press =  1.0;
            }
	  else
            {
	      rho   =  3.0;
	      u1    =  0.0;
	      u2    =  0.0;
	      u3    =  0.0;
	      press =  3.0;
            }
	}
      
      double energy = press/(gamma-1.0e0) 
	+ 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);
      
      qvals.set(i,1, rho    );
      qvals.set(i,2, rho*u1 );
      qvals.set(i,3, rho*u2 );
      qvals.set(i,4, rho*u3 );
      qvals.set(i,5, energy );
    }
}
