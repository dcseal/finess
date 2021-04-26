#include <cmath>
#include <stdexcept>
#include "assert.h"
#include "tensors.h"
#include "WenoReconstruct.h"
#include "IniParams.h"

// All purpose routine for computing a conservative WENO reconstruction that
// is based on matching polynomials with cell averages.  Upon taking
// differences of these values, you get a high-order approximation to the
// derivative.
//
// For example, g_x( x_i ) \approx ( g_{i+1/2} - g_{i-1/2} ) / dx.
//
// Input:
//
//      g( 1:meqn, 1:ws ) - list of meqn functions to be differentiated. ws =
//                          size of stencil under consideration.  ws =
//                          space_order.
//
// Output:
//
//      g_reconst( 1:meqn, 1 ) - The reconstructed value of g evaluated at 
//                          the 'right' half of the stencil, i+1/2.  
//
//      To get the other value at i-1/2, reverse the stencil, and call 
//      this same function again.
//
//
// For example, in WENO5, one passes in the following stencil:
//
//     u = { u_{i-2}, u_{i-1}, u_i, u_{i+1}, u_{i+2} },
//
// and then reconstructs the value u_{i+1/2} with this method.  To get the
// value u_{i-1/2}, pass in the following stencil:
//
//     u_{i-1/2} = { u_{i+2}, u_{i+1}, u_i, u_{i-1}, u_{i-2} }.
//
//void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst)
reconstruct_t GetWenoReconstruct()
{

    // TODO - implement WENO3 and WENO-Z version
    if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 5)
        return &WenoReconstruct_JS5;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 7)
        return &WenoReconstruct_JS7;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 9)
        return &WenoReconstruct_JS9;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 11)
        return &WenoReconstruct_JS11;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::FD && global_ini_params.get_space_order() == 5)
        return &WenoReconstruct_FD5;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::FD && global_ini_params.get_space_order() == 7)
        return &WenoReconstruct_FD7;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::FD && global_ini_params.get_space_order() == 9)
        return &WenoReconstruct_FD9;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::FD && global_ini_params.get_space_order() == 11)
        return &WenoReconstruct_FD11;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::Z  && global_ini_params.get_space_order() == 5)
        return &WenoReconstruct_Z5;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::Z  && global_ini_params.get_space_order() == 7)
        return &WenoReconstruct_Z7;
    else if(global_ini_params.get_space_order() == 1)
        return &WenoReconstructLLF;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::Z  && global_ini_params.get_space_order() == 9)
    {
//      printf("Warning: we're not sure these are the correct coefficients for WENOZ-9!\n");
        return &WenoReconstruct_Z9;
    }
    else
        throw(std::logic_error("Requested WENO Reconstruction not implemented."));
}

// ---------------------------------------------------- // 
// SECTION: WENO-JS reconstructions                     //
// ---------------------------------------------------- // 

// Fifth-order Jiang and Shu WENO reconstruction
//
//  See: G.-S. Jiang and C.-W. Shu, "Efficient Implementation of Weighted ENO 
//       Schemes". J. Comput. Phys. 126 (1996), pp. 202-228. 
//
//  or Section 2.2 of the recent SIAM review paper:
//
//       C.-W. Shu, "High Order Weighted Essentially Nonoscillatory Schemes for
//                   Convection Dominated Problems", SIAM Review, Vol. 51, 
//                   No. 1, pp 82--126, (2009).
//
// See also: WenoReconstruct.  This routine requires a stencil of lenght 5 to
// be passed in.
void WenoReconstruct_JS5( const dTensor2& g, dTensor2& g_reconst )
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


    const double eps         = global_ini_params.get_epsilon(); 
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    const int meqn = g.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // Compute smoothness indicators (identical for left/right values):
        beta0 =(13./12.)*pow(uim2-2*uim1+ui,2)+0.25*pow(uim2-4*uim1+3*ui,2);
        beta1 =(13./12.)*pow(uim1-2*ui+uip1,2)+0.25*pow(uim1-uip1,2);
        beta2 =(13./12.)*pow(ui-2*uip1+uip2,2)+0.25*pow(3*ui-4*uip1+uip2,2);
        
        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
        u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
        u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;
        
        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*pow(eps+beta0,-power_param);
        omt1 = g1*pow(eps+beta1,-power_param);
        omt2 = g2*pow(eps+beta2,-power_param);
        omts = omt0+omt1+omt2;

        // # Return 5th-order conservative reconstruction
        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3)/omts );

    }

}

// 7th-order WENO reconstruction
//
// See: D.S. Balsara and C.-W. Shu, "Monotonicity Preserving Weighted 
//      Essentially Non-oscillatory Schemes with Increasingly High Order of 
//      Accuracy". J. Comput. Phys. 160 (2000), pp. 405-452.
//
// See also: WenoReconstruct.  This routine requires a stencil of length 7 to
// be passed in.
void WenoReconstruct_JS7( const dTensor2& g, dTensor2& g_reconst )
{
    assert_eq(g.getsize(2), 7);
    double uim3, uim2, uim1, ui, uip1, uip2, uip3;
    double u1, u2, u3, u4;
    
    double beta0, beta1, beta2, beta3;
    double omt0, omt1, omt2, omt3, omts;

    const double g0 = 1./35., g1 = 12./35., g2 = 18./35., g3 = 4./35.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

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
                uim1*(11003.*uim1 -  9402.*ui  ) + 2107. * ui*ui;
        beta1 = uim2*(  267.*uim2 - 1642.*uim1   + 1602.*ui - 494.*uip1) + 
                uim1*( 2843.*uim1 - 5966.*ui     + 1922.*uip1) +
                ui*(   3443.*ui   - 2522.*uip1 ) + 547. * uip1*uip1;
        beta2 = uim1*( 547.*uim1 - 2522.*ui     + 1922.*uip1 - 494.*uip2) + 
                ui  *(3443.*ui   - 5966.*uip1   + 1602.*uip2 )   + 
                uip1*(2843.*uip1 - 1642.*uip2 ) + 267.*uip2*uip2;
        beta3 = ui*  ( 2107.*ui   -  9402.*uip1   + 7042.*uip2 - 1854.*uip3 ) + 
                uip1*(11003.*uip1 - 17246.*uip2   + 4642.*uip3 )              + 
                uip2*( 7043.*uip2 -  3882.*uip3 ) + 547.*uip3*uip3;

        u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui;
        u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1;
        u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2;
        u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3;

        omt0 = g0*pow(eps+beta0, -power_param);
        omt1 = g1*pow(eps+beta1, -power_param);
        omt2 = g2*pow(eps+beta2, -power_param);
        omt3 = g3*pow(eps+beta3, -power_param);
        omts = omt0 + omt1 + omt2 + omt3;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4)/omts);
    }
}

// 9th-order WENO reconstruction
//
// See also: WenoReconstruct.  This routine requires a stencil of length 9 to
// be passed in.
void WenoReconstruct_JS9( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 9);
    double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4;
    double u1, u2, u3, u4, u5;
    
    double beta0, beta1, beta2, beta3, beta4;
    double omt0, omt1, omt2, omt3, omt4, omts;

    const double g0 = 1.0/126.0, g1 = 10./63., g2 = 10./21., g3 = 20./63., g4 = 5./126.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++)
    {

        uim4 = g.get(m, 1);
        uim3 = g.get(m, 2);
        uim2 = g.get(m, 3);
        uim1 = g.get(m, 4);
        ui   = g.get(m, 5);
        uip1 = g.get(m, 6);
        uip2 = g.get(m, 7);
        uip3 = g.get(m, 8);
        uip4 = g.get(m, 9);

        beta0 = uim4*(22658.*uim4 - 208501.*uim3 + 364863.*uim2 - 288007.*uim1 + 86329.*ui) +
                uim3*(482963.*uim3 - 1704396.*uim2 + 1358458.*uim1 - 411487. * ui) +
                uim2*(1521393.*uim2 - 2462076.*uim1 + 758823.*ui) + 
                uim1*(1020563.*uim1 - 649501.*ui) +
                107918.*pow(ui, 2);
        beta1 = uim3*(6908.*uim3  - 60871.*uim2 + 99213.*uim1 - 70237.*ui + 18079.*uip1) +
                uim2*(138563.*uim2 - 464976.*uim1 + 337018.*ui - 88297.*uip1) +
                uim1*(406293.*uim1 - 611976*ui + 165153.*uip1) +
                ui*(242723.*ui - 140251.*uip1) + 
                22658.*pow(uip1, 2);
        beta2 = uim2*(6908.*uim2 - 51001.*uim1 + 67923.*ui - 38947.*uip1 + 8209.*uip2) +
                uim1*(104963.*uim1 - 299076.*ui + 179098.*uip1 - 38947.*uip2) +
                ui*(231153.*ui - 299076.*uip1 + 67923.*uip2) +
                uip1*(104963.*uip1 - 51001.*uip2) +
                6908.*pow(uip2, 2);
        beta3 = uim1*(22658.*uim1 - 140251.*ui + 165153.*uip1 - 88297.*uip2 + 18079.*uip3) +
                ui*(242723.*ui - 611976.*uip1 + 337018.*uip2 - 70237.*uip3) +
                uip1*(406293.*uip1 - 464976.*uip2 + 99213.*uip3) +
                uip2*(138563.*uip2 - 60871.*uip3) +
                6908.*pow(uip3, 2);
        beta4 = ui*(107918.*ui - 649501.*uip1 + 758823.*uip2 - 411487.*uip3 + 86329.*uip4) +
                uip1*(1020563.*uip1 - 2462076.*uip2 + 1358458.*uip3 - 288007*uip4) + 
                uip2*(1521393.*uip2 - 1704396.*uip3 + 364863.*uip4) +
                uip3*(482963.*uip3 - 208501.*uip4) +
                22658.*pow(uip4, 2);


        u1 = (1./5.)*uim4   + (-21./20.)*uim3 + (137./60.)*uim2 + (-163./60.)*uim1 + (137./60.)*ui;
        u2 = (-1./20.)*uim3 + (17./60.)*uim2  + (-43./60.)*uim1 + (77./60.)*ui     + (1./5.)*uip1;
        u3 = (1./30.)*uim2  + (-13./60.)*uim1 + (47./60.)*ui    + (9./20.)*uip1    + (-1./20.)*uip2;
        u4 = (-1./20.)*uim1 + (9./20.)*ui     + (47./60.)*uip1  + (-13./60.)*uip2  + (1./30.)*uip3;
        u5 = (1./5.)*ui     + (77./60.)*uip1  + (-43./60.)*uip2 + (17./60.)*uip3   + (-1./20.)*uip4;


        omt0 = g0*pow(eps+beta0, -power_param);
        omt1 = g1*pow(eps+beta1, -power_param);
        omt2 = g2*pow(eps+beta2, -power_param);
        omt3 = g3*pow(eps+beta3, -power_param);
        omt4 = g4*pow(eps+beta4, -power_param);
        omts = omt0 + omt1 + omt2 + omt3 + omt4;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4 + omt4*u5)/omts);
    }
}

// 11th-order WENO reconstruction (Jiang and Shu weights)
void WenoReconstruct_JS11( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 11);
    double uim5, uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4, uip5;
    double u1, u2, u3, u4, u5, u6;
    
    double beta0, beta1, beta2, beta3, beta4, beta5, beta6;
    double omt0, omt1, omt2, omt3, omt4, omt5, omt6, omts;

    const double g0=1.0/462., g1=5./77, g2=25./77., g3=100./231., g4=25./154., g5=1./77.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++)
    {

        uim5 = g.get(m, 1);
        uim4 = g.get(m, 2);
        uim3 = g.get(m, 3);
        uim2 = g.get(m, 4);
        uim1 = g.get(m, 5);
        ui   = g.get(m, 6);
        uip1 = g.get(m, 7);
        uip2 = g.get(m, 8);
        uip3 = g.get(m, 9);
        uip4 = g.get(m,10);
        uip5 = g.get(m,11);

        beta0 = uim5*(  1152561.*uim5 - 12950184. *uim4 + 29442256. *uim3 - 33918804. * uim2 + 19834350.*uim1 - 4712740.*ui ) 
              + uim4*( 36480687.*uim4 - 166461044.*uim3 + 192596472.*uim2 - 113206788.* uim1 + 27060170.*ui )
              + uim3*(190757572.*uim3 - 444003904.*uim2 + 262901672.*uim1 - 63394124. * ui )
              + uim2*(260445372.*uim2 - 311771244.*uim1 + 76206736. *ui )
              + uim1*(94851237. *uim1 - 47460464.*ui )  + 6150211.*ui*ui;

        beta1 = uim4*(  271779.* uim4 - 3015728.*uim3  + 6694608.  * uim2 - 7408908.  *  uim1 + 4067018.*ui - 880548.*uip1 )
              + uim3*( 8449957.* uim3 - 37913324.*uim2 + 42405032. * uim1 - 23510468. * ui    + 5134574.*uip1 )
              + uim2*(43093692.* uim2 - 97838784.*uim1 + 55053752. * ui   - 12183636. * uip1 )
              + uim1*(56662212.* uim1 - 65224244.*ui   + 14742480. * uip1 )
              + ui  *(19365967.* ui   - 9117992.*uip1 ) + 1152561. *uip1*uip1;

        beta2  = uim3*(139633.  *uim3 - 1429976.  *uim2 + 2863984. *uim1 - 2792660.*ui   + 1325006.*uip1 - 245620.*uip2 ) 
               + uim2*(3824847. *uim2 - 15880404. *uim1 + 15929912.*ui   - 7727988.*uip1 + 1458762.*uip2 ) 
               + uim1*(17195652.*uim1 - 35817664. *ui   + 17905032.*uip1 - 3462252.*uip2 ) 
               + ui  *(19510972.*ui   - 20427884. *uip1 + 4086352. *uip2 ) 
               + uip1*(5653317. *uip1 - 2380800.  *uip2 ) + 271779.*uip2*uip2;

        beta3 = 
              uim2*(  271779.* uim2 - 2380800.  *uim1 + 4086352. *ui   - 3462252.*uip1 + 1458762.*uip2 - 245620.* uip3 ) 
            + uim1*(5653317. * uim1 - 20427884. *ui    + 17905032.*uip1 - 7727988.*uip2 + 1325006.*uip3 ) 
            + ui  *(19510972.* ui   - 35817664. *uip1  + 15929912.*uip2 - 2792660.*uip3 ) 
            + uip1*(17195652.* uip1 - 15880404. *uip2  + 2863984. *uip3 ) 
            + uip2*(3824847. * uip2 - 1429976.  *uip3 ) 
            + 139633.*uip3*uip3;

        beta4 = 
            uim1*(1152561. * uim1 - 9117992. * ui   + 14742480.* uip1 - 12183636.*uip2 + 5134574.*uip3 - 880548.*uip4 ) 
          + ui  *(19365967.* ui   - 65224244.* uip1 + 55053752.* uip2 - 23510468.*uip3 + 4067018.*uip4 ) 
          + uip1*(56662212.* uip1 - 97838784.* uip2 + 42405032.* uip3 - 7408908. *uip4 ) 
          + uip2*(43093692.* uip2 - 37913324.* uip3 + 6694608. * uip4 ) 
          + uip3*(8449957.*  uip3 - 3015728. * uip4 ) 
          + 271779.*uip4*uip4;

        beta5 = 
          ui  *(6150211.*   ui   - 47460464. *uip1 + 76206736. *uip2 - 63394124. *uip3 + 27060170.*uip4 - 4712740.*uip5 ) 
        + uip1*(94851237.*  uip1 - 311771244.*uip2 + 262901672.*uip3 - 113206788.*uip4 + 19834350.*uip5 ) 
        + uip2*(260445372.* uip2 - 444003904.*uip3 + 192596472.*uip4 - 33918804. *uip5 ) 
        + uip3*(190757572.* uip3 - 166461044.*uip4 + 29442256. *uip5 ) 
        + uip4*(36480687.*  uip4 - 12950184. *uip5 ) 
        + 1152561.*uip5*uip5;

        u1 = (-1./6.)*uim5 + (31./30.)*uim4 + (-163./60.)*uim3 + (79./20.)*uim2 + (-71./20.)*uim1 + (49./20.)*ui;  
        u2 = (1./30.)*uim4 + (-13./60.)*uim3 + (37./60.)*uim2 + (-21./20.)*uim1 + (29./20.)*ui + (1./6.)*uip1;  
        u3 = (-1./60.)*uim3 + (7./60.)*uim2 + (-23./60.)*uim1 + (19./20.)*ui + (11./30.)*uip1 + (-1./30.)*uip2;  
        u4 = (1./60.)*uim2 + (-2./15.)*uim1 + (37./60.)*ui   + (37./60.)*uip1 + (-2./15.)*uip2 + (1./60.)*uip3;  
        u5 = (-1./30.)*uim1 + (11./30.)*ui  + (19./20.)*uip1 + (-23./60.)*uip2 + (7./60.)*uip3 + (-1./60.)*uip4;  
        u6 = (1./6.)*ui   + (29./20.)*uip1 + (-21./20.)*uip2 + (37./60.)*uip3 + (-13./60.)*uip4 + (1./30.)*uip5;  

        omt0 = g0*pow(eps+beta0, -power_param);
        omt1 = g1*pow(eps+beta1, -power_param);
        omt2 = g2*pow(eps+beta2, -power_param);
        omt3 = g3*pow(eps+beta3, -power_param);
        omt4 = g4*pow(eps+beta4, -power_param);
        omt5 = g5*pow(eps+beta5, -power_param);
        omts = omt0 + omt1 + omt2 + omt3 + omt4 + omt5;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4 + omt4*u5 + omt5*u6 )/omts);
    }

}

// ---------------------------------------------------- // 
// SECTION: WENO-Z reconstruction                       //
// ---------------------------------------------------- // 

// See: M. Castro, B. Costa, W.S. Don, "High order weighted essentially
//      non-oscillatory WENO-Z schemes for hyperbolic conservation laws". J. 
//      Comput.  Phys. 230 (2011), pp. 1766-1792 .
//
void WenoReconstruct_Z5( const dTensor2& g, dTensor2& g_reconst )
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

    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    const int meqn = g.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // Compute smoothness indicators (identical for left/right values):
        beta0 =(13./12.)*pow(uim2-2*uim1+ui,2)+0.25*pow(uim2-4*uim1+3*ui,2);
        beta1 =(13./12.)*pow(uim1-2*ui+uip1,2)+0.25*pow(uim1-uip1,2);
        beta2 =(13./12.)*pow(ui-2*uip1+uip2,2)+0.25*pow(3*ui-4*uip1+uip2,2);

        double tau5 = fabs( beta0 - beta2 );

        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
        u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
        u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;

        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*(1. + pow( tau5 / (eps+beta0), power_param ) );
        omt1 = g1*(1. + pow( tau5 / (eps+beta1), power_param ) );
        omt2 = g2*(1. + pow( tau5 / (eps+beta2), power_param ) );
        omts = omt0+omt1+omt2;

        // # Return 5th-order conservative reconstruction
        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3)/omts );

    }


}

void WenoReconstruct_Z7( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 7);
    double uim3, uim2, uim1, ui, uip1, uip2, uip3;
    double u1, u2, u3, u4;
    
    double beta0, beta1, beta2, beta3;
    double omt0, omt1, omt2, omt3, omts;

    const double g0 = 1./35., g1 = 12./35., g2 = 18./35., g3 = 4./35.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

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

        double tau5 = abs( beta0 + 3.*(beta1-beta2) - beta3 );

        u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui;
        u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1;
        u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2;
        u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3;

        omt0 = g0*(1. + pow( tau5 / (eps+beta0), power_param ) );
        omt1 = g1*(1. + pow( tau5 / (eps+beta1), power_param ) );
        omt2 = g2*(1. + pow( tau5 / (eps+beta2), power_param ) );
        omt3 = g3*(1. + pow( tau5 / (eps+beta3), power_param ) );
        omts = omt0 + omt1 + omt2 + omt3;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4)/omts);
    }

}

void WenoReconstruct_Z9( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 9);
    double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4;
    double u1, u2, u3, u4, u5;
    
    double beta0, beta1, beta2, beta3, beta4;
    double omt0, omt1, omt2, omt3, omt4, omts;

    const double g0 = 1.0/126.0, g1 = 10./63., g2 = 10./21., g3 = 20./63., g4 = 5./126.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++)
    {

        uim4 = g.get(m, 1);
        uim3 = g.get(m, 2);
        uim2 = g.get(m, 3);
        uim1 = g.get(m, 4);
        ui   = g.get(m, 5);
        uip1 = g.get(m, 6);
        uip2 = g.get(m, 7);
        uip3 = g.get(m, 8);
        uip4 = g.get(m, 9);

        beta0 = uim4*(22658.*uim4 - 208501.*uim3 + 364863.*uim2 - 288007.*uim1 + 86329.*ui) +
                uim3*(482963.*uim3 - 1704396.*uim2 + 1358458.*uim1 - 411487. * ui) +
                uim2*(1521393.*uim2 - 2462076.*uim1 + 758823.*ui) + 
                uim1*(1020563.*uim1 - 649501.*ui) +
                107918.*pow(ui, 2);
        beta1 = uim3*(6908.*uim3  - 60871.*uim2 + 99213.*uim1 - 70237.*ui + 18079.*uip1) +
                uim2*(138563.*uim2 - 464976.*uim1 + 337018.*ui - 88297.*uip1) +
                uim1*(406293.*uim1 - 611976*ui + 165153.*uip1) +
                ui*(242723.*ui - 140251.*uip1) + 
                22658.*pow(uip1, 2);
        beta2 = uim2*(6908.*uim2 - 51001.*uim1 + 67923.*ui - 38947.*uip1 + 8209.*uip2) +
                uim1*(104963.*uim1 - 299076.*ui + 179098.*uip1 - 38947.*uip2) +
                ui*(231153.*ui - 299076.*uip1 + 67923.*uip2) +
                uip1*(104963.*uip1 - 51001.*uip2) +
                6908.*pow(uip2, 2);
        beta3 = uim1*(22658.*uim1 - 140251.*ui + 165153.*uip1 - 88297.*uip2 + 18079.*uip3) +
                ui*(242723.*ui - 611976.*uip1 + 337018.*uip2 - 70237.*uip3) +
                uip1*(406293.*uip1 - 464976.*uip2 + 99213.*uip3) +
                uip2*(138563.*uip2 - 60871.*uip3) +
                6908.*pow(uip3, 2);
        beta4 = ui*(107918.*ui - 649501.*uip1 + 758823.*uip2 - 411487.*uip3 + 86329.*uip4) +
                uip1*(1020563.*uip1 - 2462076.*uip2 + 1358458.*uip3 - 288007*uip4) + 
                uip2*(1521393.*uip2 - 1704396.*uip3 + 364863.*uip4) +
                uip3*(482963.*uip3 - 208501.*uip4) +
                22658.*pow(uip4, 2);

        // TODO - check this!
        double tau = abs( beta4 - beta0 );

        u1 = (1./5.)*uim4   + (-21./20.)*uim3 + (137./60.)*uim2 + (-163./60.)*uim1 + (137./60.)*ui;
        u2 = (-1./20.)*uim3 + (17./60.)*uim2  + (-43./60.)*uim1 + (77./60.)*ui     + (1./5.)*uip1;
        u3 = (1./30.)*uim2  + (-13./60.)*uim1 + (47./60.)*ui    + (9./20.)*uip1    + (-1./20.)*uip2;
        u4 = (-1./20.)*uim1 + (9./20.)*ui     + (47./60.)*uip1  + (-13./60.)*uip2  + (1./30.)*uip3;
        u5 = (1./5.)*ui     + (77./60.)*uip1  + (-43./60.)*uip2 + (17./60.)*uip3   + (-1./20.)*uip4;

        omt0 = g0*(1. + pow( tau / (eps+beta0), power_param ) );
        omt1 = g1*(1. + pow( tau / (eps+beta1), power_param ) );
        omt2 = g2*(1. + pow( tau / (eps+beta2), power_param ) );
        omt3 = g3*(1. + pow( tau / (eps+beta3), power_param ) );
        omt4 = g4*(1. + pow( tau / (eps+beta4), power_param ) );
        omts = omt0 + omt1 + omt2 + omt3 + omt4;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4 + omt4*u5)/omts);
    }

}

// ---------------------------------------------------- // 
// SECTION: Upwinded reconstruction with linear weights //
// ---------------------------------------------------- // 

// Conservative reconstruction based on the *linear* weights.  Do not use for
// problems with shocks.
void WenoReconstruct_FD5( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq( g.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    // double g0, g1, g2;           
    // g0 = 0.1; g1 = 0.6; g2 = 0.3;

    const int meqn = g.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // -- central finite difference reconstruction -- //
        u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
        u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
        u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;

        // 5th-order reconstruction (linear weights)
        g_reconst.set(m, 1, 0.1*u1+0.6*u2+0.3*u3 );

    }

}

void WenoReconstruct_FD7( const dTensor2& g, dTensor2& g_reconst )
{
    assert_eq(g.getsize(2), 7);
    double uim3, uim2, uim1, ui, uip1, uip2, uip3;
    double u1, u2, u3, u4;
    
    double beta0, beta1, beta2, beta3;
    double omt0, omt1, omt2, omt3, omts;

    const double g0 = 1./35., g1 = 12./35., g2 = 18./35., g3 = 4./35.;

    const int meqn = g.getsize(1);
    
    for(int m = 1; m <= meqn; m++)
    {

        uim3 = g.get(m, 1);
        uim2 = g.get(m, 2);
        uim1 = g.get(m, 3);
        ui   = g.get(m, 4);
        uip1 = g.get(m, 5);
        uip2 = g.get(m, 6);
        uip3 = g.get(m, 7);

        u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui;
        u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1;
        u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2;
        u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3;

        g_reconst.set(m, 1, (g0*u1 + g1*u2 + g2*u3 + g3*u4) );
   }
}

// 9th-order reconstruction (linear weights)
void WenoReconstruct_FD9( const dTensor2& g, dTensor2& g_reconst )
{
    assert_eq(g.getsize(2), 9);
    double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4;
    double u1, u2, u3, u4, u5;
    
    double beta0, beta1, beta2, beta3, beta4;
    double omt0, omt1, omt2, omt3, omt4, omts;

    const double g0 = 1.0/126.0, g1 = 10./63., g2 = 10./21., g3 = 20./63., g4 = 5./126.;

    const int meqn = g.getsize(1);
    for(int m = 1; m <= meqn; m++)
    {
        uim4 = g.get(m, 1);
        uim3 = g.get(m, 2);
        uim2 = g.get(m, 3);
        uim1 = g.get(m, 4);
        ui   = g.get(m, 5);
        uip1 = g.get(m, 6);
        uip2 = g.get(m, 7);
        uip3 = g.get(m, 8);
        uip4 = g.get(m, 9);

        u1 = (1./5.)*uim4   + (-21./20.)*uim3 + (137./60.)*uim2 + (-163./60.)*uim1 + (137./60.)*ui;
        u2 = (-1./20.)*uim3 + (17./60.)*uim2  + (-43./60.)*uim1 + (77./60.)*ui     + (1./5.)*uip1;
        u3 = (1./30.)*uim2  + (-13./60.)*uim1 + (47./60.)*ui    + (9./20.)*uip1    + (-1./20.)*uip2;
        u4 = (-1./20.)*uim1 + (9./20.)*ui     + (47./60.)*uip1  + (-13./60.)*uip2  + (1./30.)*uip3;
        u5 = (1./5.)*ui     + (77./60.)*uip1  + (-43./60.)*uip2 + (17./60.)*uip3   + (-1./20.)*uip4;

        g_reconst.set(m, 1, (g0*u1 + g1*u2 + g2*u3 + g3*u4 + g4*u5) );
    }
}

// 11th-order WENO reconstruction (Jiang and Shu weights)
void WenoReconstruct_FD11( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 11);
    double uim5, uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4, uip5;
    double u1, u2, u3, u4, u5, u6;
    
    double beta0, beta1, beta2, beta3, beta4, beta5, beta6;
    double omt0, omt1, omt2, omt3, omt4, omt5, omt6, omts;

    const double g0=1.0/462., g1=5./77, g2=25./77., g3=100./231., g4=25./154., g5=1./77.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++)
    {

        uim5 = g.get(m, 1);
        uim4 = g.get(m, 2);
        uim3 = g.get(m, 3);
        uim2 = g.get(m, 4);
        uim1 = g.get(m, 5);
        ui   = g.get(m, 6);
        uip1 = g.get(m, 7);
        uip2 = g.get(m, 8);
        uip3 = g.get(m, 9);
        uip4 = g.get(m,10);
        uip5 = g.get(m,11);

        u1 = (-1./6.)*uim5 + (31./30.)*uim4 + (-163./60.)*uim3 + (79./20.)*uim2 + (-71./20.)*uim1 + (49./20.)*ui;  
        u2 = (1./30.)*uim4 + (-13./60.)*uim3 + (37./60.)*uim2 + (-21./20.)*uim1 + (29./20.)*ui + (1./6.)*uip1;  
        u3 = (-1./60.)*uim3 + (7./60.)*uim2 + (-23./60.)*uim1 + (19./20.)*ui + (11./30.)*uip1 + (-1./30.)*uip2;  
        u4 = (1./60.)*uim2 + (-2./15.)*uim1 + (37./60.)*ui   + (37./60.)*uip1 + (-2./15.)*uip2 + (1./60.)*uip3;  
        u5 = (-1./30.)*uim1 + (11./30.)*ui  + (19./20.)*uip1 + (-23./60.)*uip2 + (7./60.)*uip3 + (-1./60.)*uip4;  
        u6 = (1./6.)*ui   + (29./20.)*uip1 + (-21./20.)*uip2 + (37./60.)*uip3 + (-13./60.)*uip4 + (1./30.)*uip5;  

        g_reconst.set(m, 1, (g0*u1 + g1*u2 + g2*u3 + g3*u4 + g4*u5 + g5*u6 ) );
    }

}


void WenoReconstructLLF( const dTensor2& g, dTensor2& g_reconst )
{
    
    const int meqn = g.getsize(1);
    for(int m = 1; m <= meqn; m++)
    {
        g_reconst.set(m, 1, g.get(m,1) );
    }

}

// ------------------------------------------------- // 
// SECTION: Central Finite difference approximations //
// ------------------------------------------------- // 

// First-derivative ( using a 5 point central stencil -> fourth-order )
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

// First-derivative non-conservative, based on WENO differentation (not
// WENO-reconstruction).
void Diff1NC( double dx, const dTensor2& g, dTensor1& fx )
{

    assert_eq( g.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    double g0, g1, g2;           
    g0 = 1./6.; g1 = 2./3.; g2 = 1./6.;

    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    const int meqn = g.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // Compute smoothness indicators (identical for left/right values):
        beta0 = pow(uim2-2.*uim1+ui,2);
        beta1 = pow(uim1-2.*ui+uip1,2);
        beta2 = pow(ui-2.*uip1+uip2,2);

        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 0.5  )*uim2 - (2.   )*uim1 + ( 1.5  )*ui;
        u2 = (-0.5  )*uim1 + (0.   )*ui   + ( 0.5  )*uip1;
        u3 = (-1.5  )*ui   + (2.   )*uip1 - ( 0.5  )*uip2;
        
        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*pow(eps+beta0,-power_param);
        omt1 = g1*pow(eps+beta1,-power_param);
        omt2 = g2*pow(eps+beta2,-power_param);
        omts = omt0+omt1+omt2;

        // # Return 5th-order approximation to derivative
        fx.set(m, ( (omt0*u1 + omt1*u2 + omt2*u3)/omts) / dx );

    }

}


// Central Finite difference approximations:
//
// Second-derivative (using a 5 point central stencil)
void Diff2( double dx, const dTensor2& f, dTensor1& fxx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils
    assert_eq( f.getsize( 2 ), 5 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = ( -f.get( m, 1 ) - f.get( m, 5 ) )*(1.0/12.0);
        tmp       += (  f.get( m, 2 ) + f.get( m, 4 ) )*(4.0/ 3.0);
        tmp       += ( -f.get( m, 3 )                 )*(5.0/ 2.0);
        fxx.set( m, tmp / (dx*dx) );
    }

}

// First-derivative non-conservative, based on WENO differentation (not
// WENO-reconstruction).
void Diff2NC( double dx, const dTensor2& f, dTensor1& fxx )
{

    assert_eq( f.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    const double g0 = -1./12.; 
    const double g1 =  7./6.; 
    const double g2 = -1./12.;

    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    const int meqn = f.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = f.get(m,1);
        uim1 = f.get(m,2);
        ui   = f.get(m,3);
        uip1 = f.get(m,4);
        uip2 = f.get(m,5);

        // Compute smoothness indicators (identical for left/right values):
        beta0 = pow(uim2-2.*uim1+ui,2);
        beta1 = pow(uim1-2.*ui+uip1,2);
        beta2 = pow(ui-2.*uip1+uip2,2);

        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 1.0  )*uim2 - (2.   )*uim1 + ( 1.0  )*ui;
        u2 = ( 1.0  )*uim1 - (2.   )*ui   + ( 1.0  )*uip1;
        u3 = ( 1.0  )*ui   - (2.   )*uip1 + ( 1.0  )*uip2;
        
        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*pow(eps+beta0,-power_param);
        omt1 = g1*pow(eps+beta1,-power_param);
        omt2 = g2*pow(eps+beta2,-power_param);
        omts = omt0+omt1+omt2;

        // # Return 5th-order approximation to derivative
        fxx.set(m, ( (omt0*u1 + omt1*u2 + omt2*u3)/omts) / (dx*dx) );

    }

}


