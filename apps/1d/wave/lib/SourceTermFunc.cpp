#include "tensors.h"

// Source term function psi in the hyperbolic balance law:
//
//            q_t + f_x = psi,
//
// Input:
//
//       xpts( 1:numpts )           - The x-coordinates for a list of points
//      qvals( 1:numpts, 1:meqn )   - Solution q at each point.
//    auxvals( 1:numpts, 1:maux )   - The auxilary function at each point.
//
// Output:
//
//     psi( 1:numpts, 1:meqn )      - The source term evaluated at each point.
//
// See also: FluxFunc.
void SourceTermFunc(const dTensor1& xpts, 
                    const dTensor2& qvals, 
                    const dTensor2& auxvals,
                    dTensor2& source )
{
   
    // We are not using any source terms for now
    
    // Loop over each point
    int numpts = xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i);

        source.set(i, 1,   0.0 );

    }
}
