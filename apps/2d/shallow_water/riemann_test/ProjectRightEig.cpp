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
    
    // Project onto right eigenvectors
    for (k=1; k<=(kmax-1); k++)
    {
        Qvals.set(1, k, Wvals.get(1,k) + Wvals.get(3,k)  );

        Qvals.set(mu,k, (u1-sqrt(h))*Wvals.get(1,k) + (u1+sqrt(h))*Wvals.get(3,k) );

        Qvals.set(mv,k,  u2*Wvals.get(1,k) + Wvals.get(2,k) + u2*Wvals.get(3,k) );
    }
}
