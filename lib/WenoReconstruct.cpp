#include <cmath>
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
//      diff_g( 1:meqn, 1 ) - The derivative of g evaluated at the 'right' half
//                          of the stencil, i+1/2.  To get the same derivative
//                          at i-1/2, reverse the stencil, and call this same
//                          function again.
//
// For example, in WENO5, one passes in the following stencil:
//
//     u = { u_{i-2}, u_{i-1}, u_i, u_{i+1}, u_{i+2} },
//
// and then reconstructs the value u_{i+1/2} with this method.
void WenoReconstruct( const dTensor2& g, dTensor2& diff_g )
{

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    double g0, g1, g2;           
    g0 = 0.1; g1 = 0.6; g2 = 0.3;

    const int meqn = g.getsize(1);
    assert_eq( g.getsize(2), 5 );

// TODO - make this a user-defined input (add a [weno] section to the
// parameters file)
const double eps = 1.0e-12;  

    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // -- central finite difference reconstruction -- //
//      u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
//      u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
//      u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;

//      // 5th-order reconstruction (linear weights)
//      diff_g.set(m, 1, 0.1*u1+0.6*u2+0.3*u3 );

        // -- Fifth-order Jiang and Shu WENO reconstruction -- //

// TODO - finish filling this in 

        // Compute smoothness indicators (identical for left/right values):
        beta0 =(13./12.)*pow(uim2-2*uim1+ui,2)+0.25*pow(uim2-4*uim1+3*ui,2);
        beta1 =(13./12.)*pow(uim1-2*ui+uip1,2)+0.25*pow(uim1-uip1,2);
        beta2 =(13./12.)*pow(ui-2*uip1+uip2,2)+0.25*pow(3*ui-4*uip1+uip2,2);
        
        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
        u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
        u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;
        
        // Get linear weights and regularization parameter
        // gamma = [0.1, 0.6, 0.3]
        // eps   = cls._eps
        
        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*pow(eps+beta0,-2);
        omt1 = g1*pow(eps+beta1,-2);
        omt2 = g2*pow(eps+beta2,-2);
        omts = omt0+omt1+omt2;

        // # Return 5th-order conservative reconstruction
        // return om[0]*u1 + om[1]*u2 + om[2]*u3
        diff_g.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3)/omts );

    }

}
