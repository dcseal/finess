#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// This function is called after a full time step. That is, after all stages
// have been completed.
//
// The default behaviour is to do nothing.
void AfterFullTimeStep(double dt, const StateVars& Qold, StateVars& Qnew )
{

    dTensorBC2& qnew            = Qnew.ref_q();
    dTensorBC2& auxnew          = Qnew.ref_aux();
    const double tnew           = Qnew.get_t();

    const dTensorBC2& qold      = Qold.const_ref_q();
    const dTensorBC2& auxold    = Qold.const_ref_aux();
    const double told           = Qold.get_t();

}
