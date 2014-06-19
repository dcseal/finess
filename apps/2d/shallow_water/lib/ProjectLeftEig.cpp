#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig( int ixy, const dTensor1& Aux_ave,
        const dTensor1& Q_ave, const dTensor2& Qvals, dTensor2& Wvals)
{    

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

    // Project onto left eigenvectors
    for(int m=1; m <= mpts; m++ )
    {
        Wvals.set(1,m, (sqrt(h)+u1)/(2.0*sqrt(h))*Qvals.get(1,m) - (0.5/sqrt(h))*Qvals.get(mu,m) );

        Wvals.set(2,m, -u2*Qvals.get(1,m) + Qvals.get(mv,m) );

        Wvals.set(3,m, (sqrt(h)-u1)/(2.0*sqrt(h))*Qvals.get(1,m) + (0.5/sqrt(h))*Qvals.get(mu,m) );
    }
}
