#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig( int ixy, const dTensor1& Aux_ave,
        const dTensor1& Q_ave, const dTensor2& Qvals, dTensor2& Wvals)
{    
    int m,k;
    int meqn = Qvals.getsize(1);
    int kmax = Qvals.getsize(2)+1;
    double h,u1,u2;
    int mu,mv;

    // Direction
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
    h  = Q_ave.get(1);
    u1 = Q_ave.get(mu)/h;
    u2 = Q_ave.get(mv)/h;

    // Project onto left eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Wvals.set(1,k, (sqrt(h)+u1)/(2.0*sqrt(h))*Qvals.get(1,k) - (0.5/sqrt(h))*Qvals.get(mu,k) );

        Wvals.set(2,k, -u2*Qvals.get(1,k) + Qvals.get(mv,k) );

        Wvals.set(3,k, (sqrt(h)-u1)/(2.0*sqrt(h))*Qvals.get(1,k) + (0.5/sqrt(h))*Qvals.get(mu,k) );
    }
}
