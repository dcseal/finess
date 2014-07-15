#include "tensors.h"

// *TEMPLATE*
//
// Function that is called before a full time step
void BeforeFullTimeStep(double dt, 
        dTensorBC4& auxold, dTensorBC4& aux, 
        dTensorBC4& qold,   dTensorBC4& q)
{
    const int   mx   = q.getsize(1);
    const int   my   = q.getsize(2);
    const int   mz   = q.getsize(3);
    const int meqn   = q.getsize(4);

    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(4);
}
