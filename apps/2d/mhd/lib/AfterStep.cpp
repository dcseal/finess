#include <cassert>
#include <cmath>
#include <iostream>

#include "tensors.h"
#include "StateVars.h"
#include "IniParams.h"
#include "app_util.h"

// *TEMPLATE*
//
// Function that is called after each stage
void AfterStep(double dt, StateVars& Q )
{
    dTensorBC3&   q = Q.ref_q();
    const dTensorBC3& aux = Q.const_ref_aux();

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int    mbc = global_ini_params.get_mbc();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

    const double t = Q.get_t();

#pragma omp parallel for
    for(int i = 1; i <= mx; ++i){
        for(int j = 1; j <= my; ++j){
            const double B1 = 
                1.0 / (12.0 * dy)
                * (aux.get(i, j - 2, 1) - 8*aux.get(i, j-1, 1) 
                        + 8*aux.get(i, j+1, 1) - aux.get(i, j+2, 1));
            const double B2 =
                -1.0 / (12.0 * dx)
                * (aux.get(i-2, j, 1) - 8*aux.get(i-1, j, 1)
                        + 8*aux.get(i+1, j, 1) - aux.get(i+2, j, 1));

            double x, y;
            ij_to_xy(xlow, ylow, dx, dy, i, j, x, y);

            q.set(i, j, 6, B1);
            q.set(i, j, 7, B2);
        }
    } 
}
