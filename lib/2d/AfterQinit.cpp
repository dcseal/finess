#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// This function is called ONCE per simulation, and immediately after 
// the initial conditions are set.
//
// The default behaviour is to do nothing.
void AfterQinit( StateVars& Qnew )
{

//  To access the state variables, uncomment the following:
//
//  dTensorBC2& q       = Qnew.ref_q();
//  dTensorBC2& aux     = Qnew.ref_aux();
//  const double t      = Qnew.get_t();

}
