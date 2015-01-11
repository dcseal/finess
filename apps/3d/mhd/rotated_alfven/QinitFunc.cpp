#include <cmath>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    using std::cos;
    using std::sin;

    const int    numpts = xpts.getsize(1);
    const double gamma1 = global_ini_params.get_gamma() - 1.0;
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};
    const double t[] = {-sin(phi), cos(phi), 0.0};
    const double r[] = {-cos(phi)*sin(theta), -sin(phi)*sin(theta), cos(theta)};

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);

        double rho,u1,u2,u3,press,B1,B2,B3,energy;

        double xi = n[0]*x + n[1]*y + n[2]*z;
        double un = 0.0;
        double Bn = 1.0;
        double ut = 0.1*sin(2.0*M_PI * xi);
        double ur = 0.1*cos(2.0*M_PI * xi);
        double Bt = ut;
        double Br = ur;

        rho   =  1.0;
        u1    =  un * n[0] + ut * t[0] + ur * r[0];
        u2    =  un * n[1] + ut * t[1] + ur * r[1];
        u3    =  un * n[2] + ut * t[2] + ur * r[2];
        press =  0.1;
        B1    =  Bn * n[0] + Bt * t[0] + Br * r[0];
        B2    =  Bn * n[1] + Bt * t[1] + Br * r[1];
        B3    =  Bn * n[2] + Bt * t[2] + Br * r[2];

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
