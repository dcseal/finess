///@file lib/1d/AfterQinit.cpp
#include "tensors.h"

// *TEMPLATE*
//
// This function is called ONCE per simulation, after 
// the initial conditions are set
///@brief Hook function, called ONCE per simulation, after the initial conditions are set.
///
///- Currently, does nothing.
void AfterQinit(const dTensor2& node, dTensorBC2& aux, dTensorBC2& q)
{
    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);  
}
