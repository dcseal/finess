#include "tensors.h"

///@brief Function that is called before a full time step
///
///Currently does nothing.  Reserved for future use.
void BeforeFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q)
{
    const int   mx   = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(2);
}
