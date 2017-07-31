#include "tensors.h"
#include "IniParams.h"
#include <math.h>

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(const dTensor1& Aux_ave, 
        const dTensor1& Q_ave, 
        const dTensor2& Wvals,
        dTensor2& Qvals)
{    
    const int meqn = Qvals.getsize(1);
    const int mpts = Qvals.getsize(2);
    const double gamma = global_ini_params.get_gamma();

    double rho = Q_ave.get(1);
    double u1  = Q_ave.get(2)/rho;

    double press=pow(rho,gamma);
    double pPrime=gamma*pow(rho,gamma-1.0);

    double sqP = sqrt(pPrime);


    // Project onto right eigenvectors
    for (int m=1; m<=mpts; m++)
    {
      Qvals.set(1,m, (Wvals.get(2,m)+Wvals.get(1,m) ) );
      Qvals.set(2,m, ((u1-sqP)*Wvals.get(2,m)+(u1+sqP)*Wvals.get(1,m) ) );

    }
}
