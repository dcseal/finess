#include "tensors.h"

// *TEMPLATE*
//
// Function that is called after each stage
void AfterStep(double dt, dTensorBC3& aux, dTensorBC3& q)
{

    const int     mx = q.getsize(1);
    const int     my = q.getsize(2);
    const int   meqn = q.getsize(3);
    const int   maux = aux.getsize(3);  
}
