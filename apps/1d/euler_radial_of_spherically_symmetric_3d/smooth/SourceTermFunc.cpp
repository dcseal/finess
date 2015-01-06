#include "IniParams.h"
#include "tensors.h"

// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
        const dTensor2& qvals, 
        const dTensor2& auxvals,
        dTensor2& fvals)
{
    const double gamma = global_ini_params.get_gamma();
    int numpts=xpts.getsize();
    double x;
    double r;
    double rho, u, E, p;

    for (int i=1; i<=numpts; i++)
    {
        x = xpts.get(i);
        r = x;
        rho = qvals.get(i, 1);
        u   = qvals.get(i, 2)/rho;
        E   = qvals.get(i, 3);
        p   = (gamma-1.0)*(E - 0.5*rho*u*u);
//        fvals.set(i,1,  -2.0*rho * u / r );
//        fvals.set(i,2,  -2.0*rho*u*u / r );
//        fvals.set(i,3,  -2.0*u * (E + p) / r);

        fvals.set(i,1,  -2.0*rho * u / r );
        fvals.set(i,2,  -2.0*rho*u*u / r );
        fvals.set(i,3,  -2.0*u * (E + p) / r);
    }
}
