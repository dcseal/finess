#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor2& xpts, 
        dTensor2& auxvals)
{

    int numpts=xpts.getsize(1);
    double b,bx,by;

    // OPT = 1 is Shock Tube Problem in x Direction
    // OPT = 2 is Shock Tube Problem in y Direction
    int OPT = 2;

    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        if (OPT==1)
        {
            b  =  0.0; //0.5*exp(-pow(14.0*x-7.0,2));
            bx =  0.0; //-196.0*(2.0*x-1.0)*b;
            by =  0.0;
        }
        else
        {
            b  =  0.0; //0.5*exp(-pow(14.0*y-7.0,2));
            bx =  0.0;
            by =  0.0; //-196.0*(2.0*y-1.0)*b;
        }

        auxvals.set(i,1, b  );
        auxvals.set(i,2, bx );
        auxvals.set(i,3, by );
    }
}
