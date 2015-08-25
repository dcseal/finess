#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "IniParams.h"
#include "Components.h"     // used for defining indices into state variable
#include "gas05.h"
//#include "assert.h"

using namespace FiveMomentComponentID;      // see gas05.h (used for defining block indexes into two "Euler" systems and Maxwell

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
// Two-fluid plasma  - TODO this is hacked together to look at speed of light
//                          only
//
void SetWaveSpd(const dTensor1& xedge,
        const dTensor1& Ql,
        const dTensor1& Qr,
        const dTensor1& Auxl,
        const dTensor1& Auxr,
        double& s1,double& s2)
{ 

    // TODO - this is a complete hack!  7/9/2015 -DS
    s1 = -global_ini_params.get_cs_light()*global_ini_params.get_cc_factor();
    s2 = -s1;

    //  int n_offset = 0;   // offset for ions
    //  n_offset     = 5;   // offset for electrons
    //  n_offset     = 10;  // offset for ions

//  double GetSoundSpeed(int noffset, const dTensor1& xedge,
//          const dTensor1& Ql,   const dTensor1& Qr,
//          const dTensor1& Auxl, const dTensor1& Auxr );
//  assert_le( GetSoundSpeed( 0, xedge, Ql, Qr, Auxl, Auxr ), s2 );
//  assert_le( GetSoundSpeed( 5, xedge, Ql, Qr, Auxl, Auxr ), s2 );
}

/*
 * Routine used for debugging.  In the future, it might not hurt to run this
 * if only to verify that we do not allow sound speeds that are faster than
 * the speed of light.  -DS
 */

/*
double GetSoundSpeed(int n_offset, const dTensor1& xedge,
        const dTensor1& Ql,
        const dTensor1& Qr,
        const dTensor1& Auxl,
        const dTensor1& Auxr )
{

    const int n_rho = n_offset + _rho;      // Density
    const int n_M1  = n_offset + _M1 ;      // 1-comenent of momentum
    const int n_M2  = n_offset + _M2 ;      // 2-comenent of momentum
    const int n_M3  = n_offset + _M3 ;      // 3-comenent of momentum
    const int n_N   = n_offset + _N  ;      // Energy
 
    // Gas constant
    const double gamma  = global_ini_params.get_gamma();

    // Left states
    double const rhol    = Ql.get(n_rho);
    double const u1l     = Ql.get(n_M1)/rhol;
    double const u2l     = Ql.get(n_M2)/rhol;
    double const u3l     = Ql.get(n_M3)/rhol;
    double const energyl = Ql.get(n_N);
    double const pressl = (gamma-1.0e0)*(energyl-0.5e0*rhol*(u1l*u1l+u2l*u2l+u3l*u3l));
    // Right states
    double const rhor    = Qr.get(n_rho);
    double const u1r     = Qr.get(n_M1)/rhor;
    double const u2r     = Qr.get(n_M2)/rhor;
    double const u3r     = Qr.get(n_M3)/rhor;
    double const energyr = Qr.get(n_N);
    double const pressr  = (gamma-1.0e0)*(energyr-0.5e0*rhor*(u1r*u1r+u2r*u2r+u3r*u3r));

    // Average states
    double const rho    = 0.5e0*(rhol+rhor);
    double const u1     = 0.5e0*(u1l+u1r);
    double const u2     = 0.5e0*(u2l+u2r);
    double const u3     = 0.5e0*(u3l+u3r);
    double const press  = 0.5e0*(pressl+pressr);

printf("rhol, rhor       = %f %f \n", rhol, rhor );
printf("u1l, u1r         = %f %f \n", u1l, u1r );
printf("u2l, u2r         = %f %f \n", u2l, u2r );
printf("u3l, u3r         = %f %f \n", u3l, u3r );
printf("energyl, energyr = %f %f \n", energyl, energyr );
printf("pressl, pressr   = %f %f \n", pressl, pressr );

    return sqrt( fabs( gamma*press/rho ) );

    // Sound speeds
//  double const cl = sqrt(fabs(gamma*pressl/rhol));
//  double const cr = sqrt(fabs(gamma*pressr/rhor));
//  double const c  = sqrt(fabs(gamma*press/rho));


}
*/
