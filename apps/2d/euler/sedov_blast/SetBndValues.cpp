#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"

void SetBndValues( StateVars& Q );

// Wrappers for main Euler library
void SetBndValuesX( StateVars& Q )
{ 
//Separated boundary condition is not needed for this case
//The separated boundary is needed only if the corner values have different
//values
SetBndValues( Q ); }

void SetBndValuesY( StateVars& Q )
{ SetBndValues( Q ); }

// This is a user-supplied routine that sets the the boundary conditions
//
// For this problem, there are a total of two types of boundary conditions
// applied:
//
//      hard surface                : left and bottom boundary, and
//      zeroth-order extrapolation  : top and right boundary.
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
    // LEFT BOUNDARY (hard surface)
    // ***********************************************
    for (i=0; i>=(1-mbc); i--)
    {
        for (j=1; j<=my; j++)
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                double tmp = q.get(1-i,j,m);
                q.set(i,j,m, tmp );
            }

            // Flip the momentum
            q.set(i, j, 2, -q.get(i,j,2) );

        }
    }
    // ***********************************************

    // ***********************************************
    // RIGHT BOUNDARY (zeroth-order extrapolation)
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
    // BOTTOM BOUNDARY (hard surface)
    // ***********************************************
    for( j=0; j >= (1-mbc); j--)
    {
        for( i= 1; i <= mx; i++)
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                double tmp = q.get(i, 1-j, m);                    
                q.set(i, j, m, tmp );
            }

            // Flip the momentum
            q.set(i, j, 3, -q.get(i,j,3) );

        }
    }
    // ***********************************************


    // ***********************************************
    // TOP BOUNDARY (zeroth-order extrapolation)
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
    for(int i=1; i<=mbc; i++)
        for(int j=1; j<=mbc; j++)
        {
            for(int m=1; m<=meqn; m++)
                {     
                    q.set(1-i,1-j,m, q.get(i,j,m) );
                }
        }
    // ***********************************************


// TODO - cross terms!?
}


