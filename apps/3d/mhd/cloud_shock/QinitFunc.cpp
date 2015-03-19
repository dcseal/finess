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

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);


        double rho,u1,u2,u3,press,B1,B2,B3,energy;

        if(x < 0.05){
            rho = 3.86859;
            u1 = 11.2536;
            u2 = 0.0;
            u3 = 0.0;
            press = 167.345;
            B1 = 0.0;
            B2 = 2.1826182;
            B3 = -2.1826182;
        }
        else{
            if(sqrt(pow(x-0.25, 2) + pow(y-0.5, 2) + pow(z - 0.5, 2)) < 0.15){
                rho = 10.0;
            }
            else{
                rho = 1.0;
            }
            u1 = 0.0;
            u2 = 0.0;
            u3 = 0.0;
            press = 1.0;
            B1 = 0.0;
            B2 = 0.56418958;
            B3 = 0.56418958;
        }

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
