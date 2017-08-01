#include "tensors.h"
#include "IniParams.h"
#include <math.h>

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(const dTensor1& Aux_ave, 
        const dTensor1& Q_ave, 
        const dTensor2& Qvals,
        dTensor2& Wvals)
{    

    const int meqn = Qvals.getsize(1);
    const int mpts = Qvals.getsize(2);
    const double gamma = global_ini_params.get_gamma();

    double rho = Q_ave.get(1);
    double u1  = Q_ave.get(2)/rho;

    double press=pow(rho,gamma);
    double pPrime=gamma*pow(rho,gamma-1.0);

    double sqP = sqrt(pPrime);

    // Project onto left eigenvectors
    for (int m=1; m<=mpts; m++)
    {
      Wvals.set(1,m, 1.0/(2.0*sqP)*(Qvals.get(2,m)+(sqP-u1)*Qvals.get(1,m) ) );
      Wvals.set(2,m, 1.0/(2.0*sqP)*(-Qvals.get(2,m)+(sqP+u1)*Qvals.get(1,m) ) );
    }
}
