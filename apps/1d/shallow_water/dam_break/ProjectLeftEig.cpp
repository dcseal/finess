#include <cmath>
#include "tensors.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
// Input:
//
// Aux_ave( 1:maux )
//   Q_ave( 1:meqn )
//   Qvals( 1:meqn, 1:numpts )
//
// Output:
//
//   Wvals( 1:meqn, 1:numpts )
//
// See also: ProjectRightEig.
void ProjectLeftEig(const dTensor1& Aux_ave, 
		    const dTensor1& Q_ave,
		    const dTensor2& Qvals,
		    dTensor2& Wvals)
{    

    const int meqn = Qvals.getsize(1);
    const int npts = Qvals.getsize(2);

    // Average states
    double h = Q_ave.get(1);
    double u = Q_ave.get(2)/h;
    double sqh = sqrt(h);
    
    // Project onto left eigenvectors
    for(int n=1; n<= npts; n++ )
    {
        Wvals.set(1,n, ( (sqh+u)*Qvals.get(1,n) - Qvals.get(2,n) )/(2.0*sqh) );
        Wvals.set(2,n, ( (sqh-u)*Qvals.get(1,n) + Qvals.get(2,n) )/(2.0*sqh) );
    }
}
