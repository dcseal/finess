#include "dogdefs.h"
#include "dog_math.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
void SetWaveSpd(const dTensor1& nvec, 
        const dTensor1& xedge,
        const dTensor1& Ql, 
        const dTensor1& Qr,
        const dTensor1& Auxl, 
        const dTensor1& Auxr,
        double& s1,
        double& s2)
{   
    // Gas constant
    double const gamma   = global_ini_params.get_gamma();

    // Left states
    const double  rhol   = fabs(Ql.get(1));
    const double u1l     = Ql.get(2)/rhol;
    const double u2l     = Ql.get(3)/rhol;
    const double u3l     = Ql.get(4)/rhol;
    const double energyl = Ql.get(5);
    const double B1l     = Ql.get(6);
    const double B2l     = Ql.get(7);
    const double B3l     = Ql.get(8);
    const double Bm2l    = B1l*B1l + B2l*B2l + B3l*B3l;
    const double um2l    = u1l*u1l + u2l*u2l + u3l*u3l;
    const double pressl  = fabs(((gamma-1.0e0)*(energyl - 0.5*rhol*um2l - 0.5*Bm2l)));

    // Right states
    const double rhor    = fabs(Qr.get(1));
    const double u1r     = Qr.get(2)/rhor;
    const double u2r     = Qr.get(3)/rhor;
    const double u3r     = Qr.get(4)/rhor;
    const double energyr = Qr.get(5);
    const double B1r     = Qr.get(6);
    const double B2r     = Qr.get(7);
    const double B3r     = Qr.get(8);
    const double Bm2r    = B1r*B1r + B2r*B2r + B3r*B3r;
    const double um2r    = u1r*u1r + u2r*u2r + u3r*u3r;
    const double pressr  = fabs(((gamma-1.0e0)*(energyr - 0.5*rhor*um2r - 0.5*Bm2r)));

    // Average states
    const double rho    = 0.5e0*(rhol+rhor);
    const double u1     = 0.5e0*(u1l+u1r);
    const double u2     = 0.5e0*(u2l+u2r);
    const double u3     = 0.5e0*(u3l+u3r);
    const double press  = 0.5e0*(pressl+pressr);
    const double B1     = 0.5e0*(B1l+B1r);
    const double B2     = 0.5e0*(B2l+B2r);
    const double B3     = 0.5e0*(B3l+B3r);
    const double Bm2    = B1*B1 + B2*B2 + B3*B3;

    // Sound speed squared
    const double a2l = fabs(gamma*pressl/rhol);
    const double a2r = fabs(gamma*pressr/rhor);
    const double a2  = fabs(gamma*press/rho);

    const double unl = nvec.get(1)*u1l + nvec.get(2)*u2l + nvec.get(3)*u3l;
    const double unr = nvec.get(1)*u1r + nvec.get(2)*u2r + nvec.get(3)*u3r;
    const double un  = nvec.get(1)*u1  + nvec.get(2)*u2  + nvec.get(3)*u3;

    const double Bnl = nvec.get(1)*B1l + nvec.get(2)*B2l + nvec.get(3)*B3l;
    const double Bnr = nvec.get(1)*B1r + nvec.get(2)*B2r + nvec.get(3)*B3r;
    const double Bn  = nvec.get(1)*B1  + nvec.get(2)*B2  + nvec.get(3)*B3;

    const double dl = a2l + Bm2l/rhol;
    const double dr = a2r + Bm2r/rhor;
    const double d  = a2  + Bm2/rho;

    // Fast magnetosonic speeds
    const double nmag = sqrt(pow(nvec.get(1),2)+pow(nvec.get(2),2)+pow(nvec.get(3),2));
    const double cf_l = nmag*sqrt(0.5*( dl + sqrt(fabs(dl*dl - 4.0*a2l*Bnl*Bnl/rhol))));
    const double cf_r = nmag*sqrt(0.5*( dr + sqrt(fabs(dr*dr - 4.0*a2r*Bnr*Bnr/rhor))));
    const double cf   = nmag*sqrt(0.5*( d + sqrt(fabs(d*d - 4.0*a2*Bn*Bn/rho))));

    // Minimum speed
    s1 = Min(unl - cf_l, un - cf);

    // Maximum speed
    s2 = Max(unr + cf_r, un + cf);

}
