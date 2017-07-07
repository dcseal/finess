#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void SetBndValues( StateVars& Q )
{

    dTensorBC3& q   = Q.ref_q();
    dTensorBC3& aux = Q.ref_aux();

    int i,j,m;

    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();
    const int meqn = global_ini_params.get_meqn();
    const int mbc  = global_ini_params.get_mbc();
    const int maux = global_ini_params.get_maux();

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

    // ***********************************************
    // BOTTOM LEFT CORNER
    // ***********************************************
    for( i=0; i>=(1-mbc); i-- )
    for( j=0; j>=(1-mbc); j--)
    {
        // q values
        for (m=1; m<=meqn; m++)
        {
            double tmp = q.get(1,1,m);                    
            q.set(i,j,m, tmp );
        }
    }
    // ***********************************************

    // ***********************************************
    // BOTTOM RIGHT CORNER
    // ***********************************************
    for( i=(mx+1); i<=(mx+mbc); i++ )
    for( j=0; j>=(1-mbc); j--)
    {
        // q values
        for (m=1; m<=meqn; m++)
        {
            double tmp = q.get(mx,1,m);                    
            q.set(i,j,m, tmp );
        }
    }
    // ***********************************************


    // ***********************************************
    // TOP RIGHT CORNER
    // ***********************************************
    for( i=(mx+1); i<=(mx+mbc); i++ )
    for( j=(my+1); j<=(my+mbc); j++)
    {
        // q values
        for (m=1; m<=meqn; m++)
        {
            double tmp = q.get(mx,my,m);                    
            q.set(i,j,m, tmp );
        }
    }
    // ***********************************************

    // ***********************************************
    // TOP LEFT CORNER
    // ***********************************************
    for( i=0; i>=(1-mbc); i--)
    for( j=(my+1); j<=(my+mbc); j++ )
    {
        // q values
        for (m=1; m<=meqn; m++)
        {
            double tmp = q.get(1,my,m);                    
            q.set(i,j,m, tmp );
        }
    }
    // ***********************************************

}

// Wrappers for main Euler library
void SetBndValuesX( StateVars& Q )
{ SetBndValues( Q ); }

void SetBndValuesY( StateVars& Q )
{ SetBndValues( Q ); }
