#include "DogParamsCart1.h"
#include "assert.h"
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
// For example, in WENO5, one passes in the following stencil:
//
//     u = { u_{i-2}, u_{i-1}, u_i, u_{i+1}, u_{i+2} },
//
// and then reconstructs the value u_{i+1/2} with this method.
void WenoReconstruct( const dTensor2& g, dTensor1& diff_g )
{

    const double xlow = dogParamsCart1.get_xlow();
    const double dx   = dogParamsCart1.get_dx();

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    const int meqn = g.getsize(1);
    assert_eq( g.getsize(2), 5 );
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
        u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
        u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;

        // 5th-order reconstruction
        diff_g.set(m, 0.1*u1+0.6*u2+0.3*u3 );
    }

}
