#include "dogdefs.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
// For a scalar problem, there is nothing to do here.
//
void ProjectLeftEig(int ixy, 
    const dTensor1& Aux_ave, const dTensor1& Q_ave, 
    const dTensor2& Qvals,         dTensor2& Wvals)
{    

    const int meqn = Qvals.getsize(1);  // assert_eq( meqn, 1 )
    const int mpts = Qvals.getsize(2);

    // Project onto left eigenvectors
    for(int m=1; m<=mpts; m++)
    {
        Wvals.set(1,m, Qvals.get(1,m) );
    }

}
