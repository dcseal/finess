#include <cmath>
#include "dogdefs.h"
#include "constants.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

//      const double r2 = pow(x-0.40,2)+pow(y-0.50,2);
//      const double r  = sqrt(r2);
//      //qvals.set(i,1,sin(8.0*pi*x)*sin(8.0*pi*y));
//      qvals.set(i,1,sin(2.0*pi*x)*sin(2.0*pi*y));
//      /*if (r<0.3)
//      {  
//          qvals.set(i,1, pow( cos(5.0/3.0*pi*r) ,6) );  
//      }
//      else
//      {  qvals.set(i,1, 0.0 ); }
//      */

        if( x > 0.3 && x < 0.7 && y > 0.3 && y < 0.7 )
        {
            qvals.set(i,1,1.0);
        }
        else
        {
            qvals.set(i,1,0.0);
        }
    }
}
