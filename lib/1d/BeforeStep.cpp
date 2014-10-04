#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// Function that is called before each time step.  
//
// The default behaviour is to do nothing.
//
// See also: AfterStep.
void BeforeStep(double dt, StateVars& Q )
{

//  To access the state variables, uncomment the following:
//
//  dTensorBC2& q       = Qnew.ref_q();
//  dTensorBC2& aux     = Qnew.ref_aux();
//  const double t      = Qnew.get_t();

}
