#include "dogdefs.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
// For a scalar problem, there is nothing to do here.
//
void ProjectRightEig(int ixy, const dTensor1& Aux_ave, 
		     const dTensor1& Q_ave, const dTensor2& Wvals,
		     dTensor2& Qvals)
{    

    const int meqn = Qvals.getsize(1);  // assert_eq( meqn, 1 );
    const int mpts = Qvals.getsize(2);

    // Project onto right eigenvectors
    for(int m=1; m<=mpts; m++)
    {
        Qvals.set(1,m, Wvals.get(1,m) );
    }

}
