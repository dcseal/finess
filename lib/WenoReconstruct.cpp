#include <cmath>
#include <stdexcept>

#include "assert.h"
#include "tensors.h"

#include "WenoParams.h"
#include "DogParams.h"




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
//      TODO - this comment is not correct.  When you take differences of
//      these, ( g_{i+1/2} - g_{i-1/2} ) / dx, only THEN do you get a finite
//      difference approximation to g_x( x_i ).
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
static void WenoReconstruct_JS5( const dTensor2& g, dTensor2& diff_g )
{

    assert_eq( g.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    double g0, g1, g2;           
    g0 = 0.1; g1 = 0.6; g2 = 0.3;

    const int meqn = g.getsize(1);

// TODO - make this a user-defined input (add a [weno] section to the
// parameters file)
//const double eps = 1.0e-12;  
    const double eps = wenoParams.epsilon;  
    const double power_param = wenoParams.power_param;


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

// -- TODO finish filling in other options here -- //

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
        omt0 = g0*pow(eps+beta0,-power_param);
        omt1 = g1*pow(eps+beta1,-power_param);
        omt2 = g2*pow(eps+beta2,-power_param);
        omts = omt0+omt1+omt2;

        // # Return 5th-order conservative reconstruction
        // return om[0]*u1 + om[1]*u2 + om[2]*u3
        diff_g.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3)/omts );

    }

}

static void WenoReconstruct_JS7( const dTensor2& g, dTensor2& diff_g )
{
    assert_eq(g.getsize(2), 7);
    double uim3, uim2, uim1, ui, uip1, uip2, uip3;
    double u1, u2, u3, u4;
    
    double beta0, beta1, beta2, beta3;
    double omt0, omt1, omt2, omt3, omts;

    const double g0 = 1./35., g1 = 12./35., g2 = 18./35., g3 = 4./35.;

    const int meqn = g.getsize(1);
    
    const double eps = wenoParams.epsilon;  
    const double power_param = wenoParams.power_param;


    for(int m = 1; m <= meqn; m++){
        uim3 = g.get(m, 1);
        uim2 = g.get(m, 2);
        uim1 = g.get(m, 3);
        ui   = g.get(m, 4);
        uip1 = g.get(m, 5);
        uip2 = g.get(m, 6);
        uip3 = g.get(m, 7);

        beta0 = uim3*(  547.*uim3 -  3882.*uim2  + 4642.*uim1 - 1854.*ui) + 
                uim2*( 7043.*uim2 - 17246.*uim1  + 7042.*ui) +              
                uim1*(11003.*uim1 -  9402.*ui  ) + 2107. * pow(ui, 2);
        beta1 = uim2*(  267.*uim2 - 1642.*uim1   + 1602.*ui - 494.*uip1) + 
                uim1*( 2843.*uim1 - 5966.*ui     + 1922.*uip1) +
                ui*(   3443.*ui   - 2522.*uip1 ) + 547. * pow(uip1, 2);
        beta2 = uim1*( 547.*uim1 - 2522.*ui     + 1922.*uip1 - 494.*uip2) + 
                ui  *(3443.*ui   - 5966.*uip1   + 1602.*uip2 )   + 
                uip1*(2843.*uip1 - 1642.*uip2 ) + 267.*pow(uip2, 2);
        beta3 = ui*  ( 2107.*ui   -  9402.*uip1   + 7042.*uip2 - 1854.*uip3 ) + 
                uip1*(11003.*uip1 - 17246.*uip2   + 4642.*uip3 )              + 
                uip2*( 7043.*uip2 -  3882.*uip3 ) + 547.*pow(uip3, 2);
        u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui;
        u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1;
        u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2;
        u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3;

        omt0 = g0*pow(eps+beta0, -power_param);
        omt1 = g1*pow(eps+beta1, -power_param);
        omt2 = g2*pow(eps+beta2, -power_param);
        omt3 = g3*pow(eps+beta3, -power_param);
        omts = omt0 + omt1 + omt2 + omt3;

        diff_g.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4)/omts);
    }
}



// Central Finite difference approximations:
//
// First-derivative (using a 5 point central stencil)
void Diff1( double dx, const dTensor2& f, dTensor1& fx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils:
    assert_eq( f.getsize( 2 ), 5 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = (  f.get( m, 1 ) - f.get( m, 5 ) )*(1.0/12.0);
        tmp       += ( -f.get( m, 2 ) + f.get( m, 4 ) )*(2.0/ 3.0);
        fx.set( m, tmp / dx );
    }

}

// Central Finite difference approximations:
//
// First-derivative (using a 5 point central stencil)
//
// This is the scalar version of the above routine
double Diff1( double dx, 
    double f1, double f2, double f3, double f4, double f5 )
{

    double tmp = (  f1 - f5 )*(1.0/12.0);
    tmp       += ( -f2 + f4 )*(2.0/ 3.0);
    return tmp/dx;

}
// Central Finite difference approximations:
//
// Second-derivative (using a 5 point central stencil)
void Diff2( double dx, const dTensor2& f, dTensor1& fxx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils:
    assert_eq( f.getsize( 2 ), 5 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = ( -f.get( m, 1 ) - f.get( m, 5 ) )*(1.0/12.0);
        tmp       += (  f.get( m, 2 ) + f.get( m, 4 ) )*(4.0/ 3.0);
        tmp       += ( -f.get( m, 3 )                 )*(5.0/ 2.0);
        fxx.set( m, tmp / (dx*dx) );
    }

}

void WenoReconstruct(const dTensor2& g, dTensor2& diff_g){
    if(wenoParams.weno_version == WENOParams::JS && dogParams.get_space_order() == 5)
        WenoReconstruct_JS5(g, diff_g);
    else if(wenoParams.weno_version == WENOParams::JS && dogParams.get_space_order() == 7)
        WenoReconstruct_JS7(g, diff_g);
    else
        throw(std::logic_error("Requested WENO Reconstruction not implemented."));
}

