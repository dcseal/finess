#include "tensors.h"

// *TEMPLATE*
//
// This function is called ONCE per simulation, after 
// the initial conditions are set
void AfterQinit(dTensorBC2& aux, dTensorBC2& q)
{
    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);  
}
