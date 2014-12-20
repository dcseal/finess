#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// This function is called ONCE per simulation, after 
// the initial conditions are set
void AfterQinit( StateVars& Qnew )
{

    dTensorBC4& q       = Qnew.ref_q();
    dTensorBC4& aux     = Qnew.ref_aux();
    const double t      = Qnew.get_t();

}
