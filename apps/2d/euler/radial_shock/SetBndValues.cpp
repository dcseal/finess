#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void SetBndValues( dTensorBC3& aux, dTensorBC3& q )
{

    int i,j,m;

    const int mx   = dogParamsCart2.get_mx();
    const int my   = dogParamsCart2.get_my();
    const int meqn = dogParams.get_meqn();
    const int mbc  = dogParamsCart2.get_mbc();
    const int maux = dogParams.get_maux();

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------

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
                double tmp = q.get(1,j,m);
                q.set(i,j,m, tmp );
            }
            // aux values
            for (m=1; m<=maux; m++)
            {
                double tmp = q.get(1,j,m);
                aux.set(i,j,m, tmp );
            }

        }
    }
    // ***********************************************



    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for( i=(mx+1); i<=(mx+mbc); i++ )
    {
        for( j=1; j<=my; j++ )
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                double tmp = q.get(mx,j,m);                    
                q.set(i,j,m, tmp );
            }
            // aux values
            for (m=1; m<=maux; m++)
            {
                double tmp = q.get(mx,j,m);                    
                aux.set(i,j,m, tmp );
            }

        }
    }
    // ***********************************************



    // ***********************************************
    // BOTTOM BOUNDARY
    // ***********************************************
    for( j=0; j >= (1-mbc); j--)
    {
        for( i= 1; i <= mx; i++)
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                double tmp = q.get(i, 1, m);                    
                q.set(i, j, m, tmp );
            }
        }
    }
    // ***********************************************


    // ***********************************************
    // TOP BOUNDARY
    // ***********************************************
    for( j=(my+1); j<=(my+mbc); j++)
    {
        for( i= 1; i<= mx; i++ )
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                double tmp = q.get(i,my,m);                    
                q.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************

}
