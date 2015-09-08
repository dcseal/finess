#include <fstream>
#include "tensors.h"
using namespace std;

// This is a user-required routine that defines the auxilary arrays when the
// parameter maux > 0.  
//
// Euler equatons does not use this routine.  (Formerly it was used to save
// the gamma gas constant, but that has since been deprecated).
//
//
// Input:
//
//    xpts( 1:numpts )             - The x-coordinates for a list of points
//
// Output:
//
//    auxvals( 1:numpts, 1:maux )  - The vector containg auxilary values.
//
// See also: QinitFunc.
void AuxFunc(const dTensor1& xpts, 
	     dTensor2& auxvals)
{

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
    }

}
