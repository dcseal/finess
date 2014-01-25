#include "tensors.h"

// *TEMPLATE*
//
// Function that is called after a full time step (i.e., after all stages are complete)
void AfterFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q)
{
    const int   mx   = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(4);
}
