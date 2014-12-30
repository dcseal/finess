#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor2& xpts, dTensor2& auxvals)
{

    const int numpts=xpts.getsize(1);

    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i,1);
        double y = xpts.get(i,2);
        double z = xpts.get(i,3);


    }
}
