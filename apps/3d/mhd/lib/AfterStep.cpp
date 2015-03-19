#include <cassert>
#include <cmath>
#include <iostream>

#include "tensors.h"
#include "StateVars.h"
#include "IniParams.h"

inline void ijk_to_xyz(double xlow, double ylow, double zlow, 
        double dx, double dy, double dz, 
        int i, int j, int k, 
        double& x, double& y, double& z){
    x = xlow + (i - 0.5) * dx;
    y = ylow + (j - 0.5) * dy;
    z = zlow + (k - 0.5) * dz;
}


// Function that is called after each stage
void AfterStep(double dt, StateVars& Q )
{
    dTensorBC4&   q = Q.ref_q();
    const dTensorBC4& aux = Q.const_ref_aux();

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int     mz = global_ini_params.get_mz();
    const int    mbc = global_ini_params.get_mbc();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();
    const double zlow = global_ini_params.get_zlow();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();
    const double dz = global_ini_params.get_dz();

    const double t = Q.get_t();

#pragma omp parallel for
    for(int i = 1; i <= mx; ++i){
        for(int j = 1; j <= my; ++j){
            for(int k = 1; k <= mz; ++k){
                if(global_ini_params.get_constrained_transport()){
                    const double old_E = q.get(i, j, k, 5);
                    const double old_B1 = q.get(i, j, k, 6);
                    const double old_B2 = q.get(i, j, k, 7);
                    const double old_B3 = q.get(i, j, k, 8);

                    const double B1 = 
                        1.0 / (12.0 * dy)
                        *(aux.get(i, j - 2, k, 3) - 8*aux.get(i, j-1, k, 3) 
                                + 8*aux.get(i, j+1, k, 3) - aux.get(i, j+2, k, 3))
                        -1.0 / (12.0 * dz)
                        *(aux.get(i, j, k-2, 2) - 8*aux.get(i, j, k-1, 2)
                                + 8*aux.get(i, j, k+1, 2) - aux.get(i, j, k+2, 2));
                    const double B2 =
                        1.0 / (12.0 * dz)
                        *(aux.get(i, j, k-2, 1) - 8*aux.get(i, j, k-1, 1)
                                + 8*aux.get(i, j, k+1, 1) - aux.get(i, j, k+2, 1))
                        -1.0 / (12.0 * dx)
                        * (aux.get(i-2, j, k, 3) - 8*aux.get(i-1, j, k, 3)
                                + 8*aux.get(i+1, j, k, 3) - aux.get(i+2, j, k, 3));
                    const double B3 =
                        1.0 / (12.0*dx)
                        *(aux.get(i-2, j, k, 2) - 8*aux.get(i-1, j, k, 2)
                                + 8*aux.get(i+1, j, k, 2) - aux.get(i+2, j, k, 2))
                        -1.0 / (12.0*dy)
                        *(aux.get(i, j-2, k, 1) - 8*aux.get(i, j-1, k, 1)
                                + 8*aux.get(i, j+1, k, 1) - aux.get(i, j+2, k, 1));
                    double x, y, z;
                    ijk_to_xyz(xlow, ylow, zlow, dx, dy, dz, i, j, k, x, y, z);


                    q.set(i, j, k, 6, B1);
                    q.set(i, j, k, 7, B2);
                    q.set(i, j, k, 8, B3);
                }

                //                if(global_ini_params.get_mpp_limiter())
                //                    q.set(i, j, 5, old_E + 0.5*((B1*B1+B2*B2+B3*B3) - (old_B1*old_B1+old_B2*old_B2+old_B3*old_B3)));
            }
        }
    }
}
