#include <cmath>

#include "IniParams.h"

#include "tensors.h"
#include "StateVars.h"

#include "app_util.h"





// This is a user-supplied routine that sets the the boundary conditions
//
// Double periodic
//
void SetBndValues( StateVars& Q )
{
    using std::cos;
    using std::sin;

    dTensorBC3& q     = Q.ref_q();
    dTensorBC3& aux   = Q.ref_aux();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    double xlow = global_ini_params.get_xlow();
    double ylow = global_ini_params.get_ylow();
    double dx = global_ini_params.get_dx();
    double dy = global_ini_params.get_dy();

    double t = Q.get_t();

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for(int i=0; i>=(1-mbc); i--)
    {
        for(int j=1; j<=my; j++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(mx + i,j,m);
                q.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************



    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for(int i=(mx+1); i<=(mx+mbc); i++)
    {
        for(int j=1; j<=my; j++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(i - mx,j,m);
                q.set(i,j,m, tmp );
            }
           
        }
    }
    // ***********************************************



    // ***********************************************
    // BOTTOM BOUNDARY
    // ***********************************************
    for(int j=0; j>=(1-mbc); j--)
    {
        for(int i=1; i<=mx; i++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(i,my + j,m);
                q.set(i,j,m, tmp );
            }               
        }
    }
    // ***********************************************



    // ***********************************************
    // TOP BOUNDARY
    // ***********************************************
    for(int j=(my+1); j<=(my+mbc); j++)
    {
        for(int i=1; i<=mx; i++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(i,j - my,m);
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
                q.set(1-i,1-j,m, q.get(mx + 1 - i, my + 1 - j,m) );
            }
        }
    // ***********************************************


    // ***********************************************
    // BOTTOM RIGHT CORNER
    // ***********************************************
    for(int i=1; i<=mbc; i++)
        for(int j=1; j<=mbc; j++)
        {
            for(int m=1; m<=meqn; m++)
            {     
                q.set(mx+i,1-j,m, q.get(i, my + 1 - j, m) );
            }

        }
    // ***********************************************


    // ***********************************************
    // TOP RIGHT CORNER
    // ***********************************************
    for(int i=1; i<=mbc; i++)
        for(int j=1; j<=mbc; j++)
        {
            for(int m=1; m<=meqn; m++)
            {     
                q.set(mx+i,my+j,m, q.get(i, j, m) );
            }
        }
    // ***********************************************


    // ***********************************************
    // TOP LEFT CORNER
    // ***********************************************
    for(int i=1; i<=mbc; i++)
        for(int j=1; j<=mbc; j++)
        {
            for(int m=1; m<=meqn; m++)
            {     
                q.set(1-i,my+j,m, q.get(mx + 1 - i, j, m) );
            }
        }
    // ***********************************************


}
