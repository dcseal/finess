#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// This function is called ONCE per simulation, after 
// the initial conditions are set
void AfterQinit( StateVars& Qstate )
{

    dTensorBC2& q   = Qstate.ref_q();
    dTensorBC2& aux = Qstate.ref_aux();
    const double t        = Qstate.get_t();

}
