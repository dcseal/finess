#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig( int ixy, const dTensor1& Aux_ave,
        const dTensor1& Q_ave, const dTensor2& Wvals, dTensor2& Qvals)
{    

    int m,k;
    const int meqn = Qvals.getsize(1);
    const int mpts = Qvals.getsize(2);

    // Direction
    int mu,mv;
    if (ixy==1)
    {  
        mu = 2;
        mv = 3;
    }
    else
    {
        mu = 3;
        mv = 2;
    }

    // Average states
    const double h  = Q_ave.get(1);
    const double u1 = Q_ave.get(mu)/h;
    const double u2 = Q_ave.get(mv)/h;

    // Project onto right eigenvectors
    for(int m=1; m<=mpts; m++)
    {
        Qvals.set(1, m, Wvals.get(1,m) + Wvals.get(3,m)  );

        Qvals.set(mu,m, (u1-sqrt(h))*Wvals.get(1,m) + (u1+sqrt(h))*Wvals.get(3,m) );

        Qvals.set(mv,m,  u2*Wvals.get(1,m) + Wvals.get(2,m) + u2*Wvals.get(3,m) );
    }

}
