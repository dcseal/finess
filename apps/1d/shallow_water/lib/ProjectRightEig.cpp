#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
// Input:
//
// Aux_ave( 1:maux )
//   Q_ave( 1:meqn )
//   Wvals( 1:meqn, 1:numpts )
//
// Output:
//
//   Qvals( 1:meqn, 1:numpts )
//
// See also: ProjectLeftEig
void ProjectRightEig(const dTensor1& Aux_ave, 
        const dTensor1& Q_ave, 
        const dTensor2& Wvals,
        dTensor2& Qvals)
{    

    const int meqn = Qvals.getsize(1);
    const int npts = Qvals.getsize(2);

    // Average states
    const double h   = Q_ave.get(1);
    const double u   = Q_ave.get(2)/h;
    const double sqh = sqrt(h);

    // Project onto right eigenvectors
    for(int n=1; n<= npts; n++ )
    {
        Qvals.set(1,n, Wvals.get(1,n) + Wvals.get(2,n) );
        Qvals.set(2,n, (u-sqh)*Wvals.get(1,n) + (u+sqh)*Wvals.get(2,n) );
    }

}
