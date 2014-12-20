#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// Function that is called before a full time step
void BeforeFullTimeStep(double dt, StateVars& Qold, StateVars& Q)
{
    dTensorBC4& q    = Q.ref_q();
    dTensorBC4& aux  = Q.ref_aux();
    double time      = Q.get_t();

    const int   mx   = q.getsize(1);
    const int   my   = q.getsize(2);
    const int   mz   = q.getsize(3);
    const int meqn   = q.getsize(4);

    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(4);

}
