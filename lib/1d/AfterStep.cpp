#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// This function that is called after each time step
//
// The default behaviour is to do nothing.
//
// See also: BeforeStep.
void AfterStep(double dt, StateVars& Q )
{

//  To access the state variables, uncomment the following:
//
//  dTensorBC2& q       = Qnew.ref_q();
//  dTensorBC2& aux     = Qnew.ref_aux();
//  const double t      = Qnew.get_t();


}
