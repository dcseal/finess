#include "tensors.h"

// Zeroth order extrapolation boundary conditions
void SetBndValues(
        dTensorBC2& aux, 
        dTensorBC2& q)
{

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for (int i=0; i>=(1-mbc); i--)
    {        
        // q values (mirror image for all variables except momentum)
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(1-i,m);
            q.set(i, m, tmp );
        }

        // Flip the momentum
        q.set(i, 2, -q.get(i,2) );

        // aux values
        for (int m=1; m<=maux; m++)
        {
            double tmp = aux.get(1,m);

            aux.set(i,m, tmp );
        }
    }
    // ***********************************************  


    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for (int i=(mx+1); i<=(mx+mbc); i++)
    {        
        // q values
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(2*mx+1-i,m);
            q.set(i, m, tmp );
        }

        // Flip the momentum
        q.set(i, 2, -q.get(i,2) );

        // aux values
        for (int m=1; m<=maux; m++)
        {
            double tmp = aux.get(mx,m);

            aux.set(i,m, tmp );
        }
    }
    // ***********************************************

}
