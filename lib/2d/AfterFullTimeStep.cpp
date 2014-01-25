#include "tensors.h"

// *TEMPLATE*
//
// Function that is called after a full time step (i.e., after all stages are complete)
void AfterFullTimeStep(double dt,
		       dTensorBC3& auxold, dTensorBC3& aux, 
		       dTensorBC3& qold,   dTensorBC3& q)
{
    const int   mx   = q.getsize(1);
    const int   my   = q.getsize(2);
    const int meqn   = q.getsize(3);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(3);
}
