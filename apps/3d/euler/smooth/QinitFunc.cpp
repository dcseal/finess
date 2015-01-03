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

    const double gamma = global_ini_params.get_gamma();
    const int    numpts = xpts.getsize(1);

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);

        double rad = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        double rho,press,u1,u2,u3,energy;

	rho = 1.0 + 0.1*exp(-30.0*pow(rad-1.0, 2));
   	u1 = 0.0;
	u2 = 0.0;
	u3 = 0.0;
	energy = rho;

        qvals.set(i,1, rho );       // density
        qvals.set(i,2, rho*u1 );    // 1-momentum
        qvals.set(i,3, rho*u2 );    // 2-momentum
        qvals.set(i,4, rho*u3 );    // 3-momentum
        qvals.set(i,5, energy );    // energy
    }

}
