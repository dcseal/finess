#include "tensors.h"
#include "StateVars.h"

// We do not use this in the implicit solver
// Left in so that FINESS compilation/execution does not break
void SetBndValues( StateVars& Q )
{

    dTensorBC2&  q  = Q.ref_q  ();
    dTensorBC2& aux = Q.ref_aux();
    double t        = Q.get_t  ();


}
