#include "tensors.h"

// *TEMPLATE*
//
// This function is called ONCE per simulation, after 
// the initial conditions are set
void AfterQinit( dTensorBC4& aux, dTensorBC4& q )
{
    const int     mx  = q.getsize(1);
    const int     my  = q.getsize(2);
    const int     mz  = q.getsize(3);
    const int    meqn = q.getsize(4);
    const int   maux  = aux.getsize(4);  
}
