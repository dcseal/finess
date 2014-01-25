#include "tensors.h"

// Function that is called after each time step
void AfterStep(double dt, const dTensor2& node, dTensorBC2& aux, dTensorBC2& q)
{
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);    
}

// Function that is called after each time step
void AfterStep(double dt, const dTensor2& node, dTensorBC2& aux, dTensorBC2& q, void* data)
{
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);    
}
