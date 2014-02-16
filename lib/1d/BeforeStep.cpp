#include "tensors.h"

///@brief Function that is called before each time step
///
///Currently does nothing.  Reserved for future use.
void BeforeStep(double dt, const dTensor2& node, dTensorBC2& aux, dTensorBC2& q)
{
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);  
}

///@brief Function that is called before each time step
///
///Currently does nothing.  Reserved for future use.
void BeforeStep(double dt, const dTensor2& node, dTensorBC2& aux, dTensorBC2& q, 
        void* data)
{
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);  
}
