#include "dogdefs.h"
#include "tensors.h"
#include "DogParamsCart2.h"

// Function that is called after a full time step (i.e., after all stages are complete)
//
// This routine exists because the domain for this problem is weird.  Here, we
// reset the values of the function to be constant inside the wedge.  Next
// time ConstructL is called, these values get overwritten.
void AfterFullTimeStep(double dt,
		       dTensorBC3& auxold, dTensorBC3& aux, 
		       dTensorBC3& qold,   dTensorBC3& q)
{
    const int   mx   = q.getsize(1);
    const int   my   = q.getsize(2);
    const int meqn   = q.getsize(3);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(3);

    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    // Reset the domain inside the wedge to be a constant
    //
    // This will be overwritten next time ConstructL is called.
    for( int i = 1-mbc; i <= istep-1; i++ )
    for( int j = 1-mbc; j <= jstep-1; j++ )
    {
        for (int m=1; m<=meqn; m++)
        {

            q.set(i,j,m, 0.0 );
        }
    }

}
