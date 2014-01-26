#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);
    double h,u1,u2,b;

    // OPT = 1 is Shock Tube Problem in x Direction
    // OPT = 2 is Shock Tube Problem in y Direction
    const int OPT = 1;

    // Loop over grid points
    for(int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        if (OPT==1)
        {
            b  =  0.0; //0.5*exp(-pow(14.0*x-7.0,2));

            if(x>0.5)
            {
                h   =  1.0-b;
                u1  =  0.0;
                u2  =  0.0;
            }
            else
            {
                h   =  3.0-b;
                u1  =  0.0;
                u2  =  0.0;
            } 
        }
        else
        {
            //b  =  0.5*exp(-pow(14.0*y-7.0,2));
            b  =  0.0;

            if(y>0.5)
            {
                h   =  1.0-b;
                u1  =  0.0;
                u2  =  0.0;
            }
            else
            {
                h   =  3.0-b;
                u1  =  0.0;
                u2  =  0.0;
            } 
        }

        qvals.set(i, 1, h    );
        qvals.set(i, 2, h*u1 );
        qvals.set(i, 3, h*u2 );
    }
}
