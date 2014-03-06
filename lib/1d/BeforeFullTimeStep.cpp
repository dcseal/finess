#include "tensors.h"

// Function that is called before a full time step
void BeforeFullTimeStep(double dt, 
		       dTensorBC2& auxold, dTensorBC2& aux, 
		       dTensorBC2& qold,   dTensorBC2& q)
{
    const int   mx   = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(2);
}
