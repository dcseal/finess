#include <cmath>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    using std::exp;
    using std::pow;
    using std::sqrt;

    const int    numpts = xpts.getsize(1);
    const double gamma1 = global_ini_params.get_gamma() - 1.0;

    const double rhol = global_ini_params.get_rhol();
    const double unl  = global_ini_params.get_unl();
    const double utl  = global_ini_params.get_utl();
    const double url  = global_ini_params.get_url();
    const double pl   = global_ini_params.get_pl();
    const double Bnl  = global_ini_params.get_Bnl();
    const double Btl  = global_ini_params.get_Btl();
    const double Brl  = global_ini_params.get_Brl();

    const double rhor = global_ini_params.get_rhor();
    const double unr  = global_ini_params.get_unr();
    const double utr  = global_ini_params.get_utr();
    const double urr  = global_ini_params.get_urr();
    const double pr   = global_ini_params.get_pr();
    const double Bnr  = global_ini_params.get_Bnr();
    const double Btr  = global_ini_params.get_Btr();
    const double Brr  = global_ini_params.get_Brr();

    const double n[] = {0, 1, 0};
    const double t[] = {0, 0, 1};
    const double r[] = {1, 0, 0};

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);


        double rho,u1,u2,u3,press,B1,B2,B3,energy;

        double xi = y;


        double un, ut, ur, Bn, Bt, Br;
        if(xi < 0.0){
            rho = rhol;
            un = unl;
            ut = utl;
            ur = url;
            press = pl;
            Bn = Bnl;
            Bt = Btl;
            Br = Brl;
        }
        else{
            rho = rhor;
            un = unr;
            ut = utr;
            ur = urr;
            press = pr;
            Bn = Bnr;
            Bt = Btr;
            Br = Brr;
        }

        u1 = n[0]*un + t[0]*ut + r[0]*ur;
        u2 = n[1]*un + t[1]*ut + r[1]*ur;
        u3 = n[2]*un + t[2]*ut + r[2]*ur;
        B1 = n[0]*Bn + t[0]*Bt + r[0]*Br;
        B2 = n[1]*Bn + t[1]*Bt + r[1]*Br;
        B3 = n[2]*Bn + t[2]*Bt + r[2]*Br;

        energy = press/gamma1 
            + 0.5*rho*(u1*u1 + u2*u2 + u3*u3)
            + 0.5*(B1*B1 + B2*B2 + B3*B3);

        qvals.set(i,1, rho );     // density
        qvals.set(i,2, rho*u1);  // 1-momentum
        qvals.set(i,3, rho*u2 );  // 2-momentum
        qvals.set(i,4, rho*u3 );  // 3-momentum
        qvals.set(i,5, energy );  // energy
        qvals.set(i,6, B1 );      // B1
        qvals.set(i,7, B2 );      // B2
        qvals.set(i,8, B3 );      // B3      
    }

}
