#include "tensors.h"

// Function that is called before each time step
void BeforeStep(double dt, dTensorBC2& aux, dTensorBC2& q)
{
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);  
}

void BeforeStep(double dt, dTensorBC2& aux, dTensorBC2& q, 
        void* data)
{
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);  
}
