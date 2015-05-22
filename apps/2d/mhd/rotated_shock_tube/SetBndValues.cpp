#include <cmath>

#include "IniParams.h"

#include "tensors.h"
#include "StateVars.h"

#include "app_util.h"



// This is a user-supplied routine that sets the the boundary conditions
//
//
void SetBndValues( StateVars& Q )
{
    using std::cos;
    using std::sin;
    using std::tan;

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
    double angle = global_ini_params.get_angle();

    double t = Q.get_t();

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for(int i=0; i>=(1-mbc); i--)
    {
        for(int j=1-mbc; j<=my+mbc; j++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(1,j,m);
                q.set(i,j,m, tmp );
            }
            aux.set(i, j, 1,
                    2.0 * aux.get(i+1, j, 1) - aux.get(i+2, j, 1));
        }
    }
    // ***********************************************



    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for(int i=(mx+1); i<=(mx+mbc); i++)
    {
        for(int j=1-mbc; j<=my+mbc; j++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(mx, j, m);
                q.set(i,j,m, tmp );
            }
           
            aux.set(i, j, 1,
                    2.0 * aux.get(i - 1, j, 1) - aux.get(i - 2, j, 1));
        }
    }
    // ***********************************************


//    // ***********************************************
//    // BOTTOM LEFT CORNER
//    // ***********************************************
//    for(int i=1; i<=mbc; i++)
//        for(int j=1; j<=mbc; j++)
//        {
//            for(int m=1; m<=meqn; m++)
//            {     
//                q.set(1-i,1-j,m, q.get(1, 1, m) );
//            }
//            aux.set(1-i, 1-j, 1,
//                    2.0 * aux.get(1, 1-j+1, 1) - aux.get(1, 1-j+2, 1));
//        }
//    // ***********************************************
//
//
//    // ***********************************************
//    // BOTTOM RIGHT CORNER
//    // ***********************************************
//    for(int i=1; i<=mbc; i++)
//        for(int j=1; j<=mbc; j++)
//        {
//            for(int m=1; m<=meqn; m++)
//            {     
//                q.set(mx+i,1-j,m, q.get(mx, 1, m) );
//            }
//
//            aux.set(mx + i, 1-j, 1,
//                    2.0 * aux.get(mx, 1-j+1, 1) - aux.get(mx, 1-j+2, 1));
//        }
//    // ***********************************************
//
//
//    // ***********************************************
//    // TOP RIGHT CORNER
//    // ***********************************************
//    for(int i=1; i<=mbc; i++)
//        for(int j=1; j<=mbc; j++)
//        {
//            for(int m=1; m<=meqn; m++)
//            {     
//                q.set(mx+i,my+j,m, q.get(mx, my, m) );
//            }
//            aux.set(mx + i, my + j, 1,
//                    2.0 * aux.get(mx, my + j - 1, 1) - aux.get(mx, my + j - 2, 1));
//        }
//    // ***********************************************
//
//
//    // ***********************************************
//    // TOP LEFT CORNER
//    // ***********************************************
//    for(int i=1; i<=mbc; i++)
//        for(int j=1; j<=mbc; j++)
//        {
//            for(int m=1; m<=meqn; m++)
//            {     
//                q.set(1-i,my+j,m, q.get(1, my, m) );
//            }
//            aux.set(1 - i, my + j, 1,
//                    2.0 * aux.get(1, my + j - 1, 1) - aux.get(1, my + j - 2, 1));
//        }
//    // ***********************************************



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
                double tmp = q.get(i-1, j+2, m);
                q.set(i,j,m, tmp );
            }               
            aux.set(i, j, 1,
                    2.0 * aux.get(i-1, j + 2, 1) - aux.get(i-2, j + 4, 1));
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
                double tmp = q.get(i+1, j - 2, m);
                q.set(i,j,m, tmp );
            }

            aux.set(i, j, 1,
                    2.0 * aux.get(i+1, j-2, 1) - aux.get(i+2, j-4, 1));
        }
    }
    // ***********************************************


}

void SetBndValuesX( StateVars& Q )
{ SetBndValues( Q ); }

void SetBndValuesY( StateVars& Q )
{ SetBndValues( Q ); }
