#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(const dTensor1& Aux_ave, 
		    const dTensor1& Q_ave,
		    const dTensor2& Qvals,
		    dTensor2& Wvals)
{    
    int m,k;
    int meqn = Qvals.getsize(1);
    int kmax = Qvals.getsize(2)+1;
    double h,u,sqh;

    // Average states
    h = Q_ave.get(1);
    u = Q_ave.get(2)/h;
    sqh = sqrt(h);
    
    // Project onto left eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Wvals.set(1,k, ((sqh+u)*Qvals.get(1,k) - Qvals.get(2,k))/(2.0*sqh) );

        Wvals.set(2,k, ((sqh-u)*Qvals.get(1,k) + Qvals.get(2,k))/(2.0*sqh) );
    }
}
