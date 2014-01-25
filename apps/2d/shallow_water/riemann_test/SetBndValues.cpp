#include "dogdefs.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void SetBndValues(dTensorBC3& q, dTensorBC3& aux)
{

    int i,j,m;
    double tmp;
    int mx   = q.getsize(1);
    int my   = q.getsize(2);
    int meqn = q.getsize(3);
    int mbc  = q.getmbc();
    int maux = aux.getsize(3);

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------

    {

        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (i=0; i>=(1-mbc); i--)
        {
            for (j=1; j<=my; j++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(1,j,m);
                    q.set(i,j,m, tmp );
                }
            }
        }
        // ***********************************************



        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (i=(mx+1); i<=(mx+mbc); i++)
        {
            for (j=1; j<=my; j++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(mx,j,m);                    
                    q.set(i,j,m, tmp );
                }
            }
        }
        // ***********************************************



        // ***********************************************
        // BOTTOM BOUNDARY
        // ***********************************************
        for (j=0; j>=(1-mbc); j--)
        {
            for (i=(2-mbc); i<=(mx+mbc-1); i++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i,1,m);                    
                    q.set(i,j,m, tmp );
                }
            }
        }
        // ***********************************************



        // ***********************************************
        // TOP BOUNDARY
        // ***********************************************
        for (j=(my+1); j<=(my+mbc); j++)
        {
            for (i=(2-mbc); i<=(mx+mbc-1); i++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i,my,m);                    
                    q.set(i,j,m, tmp );
                }
            }
        }
        // ***********************************************

    }

}
