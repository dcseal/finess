#include "DogParamsCart1.h"
#include "tensors.h"

// All purpose routine for computing a conservative finite difference
// approximation to the derivative of the function.
//
// Input:
//
//      g( 1:meqn, 1:ws ) - list of meqn functions to be differentiated. ws =
//                          size of stencil under consideration.  ws =
//                          space_order.
//
// Output:
//
//      diff_g( 1:meqn  ) - The derivative of g evaluated at the 'right' half
//                          of the stencil, i+1/2.  To get the same derivative
//                          at i-1/2, reverse the stencil, and call this same
//                          function again.
//
void WenoReconstruct( const dTensor2& g, dTensor1& diff_g )
{

    const double xlow = dogParamsCart1.get_xlow();
    const double dx   = dogParamsCart1.get_dx();

}
