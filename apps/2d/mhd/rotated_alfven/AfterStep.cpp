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

    const double angle = global_ini_params.get_angle();
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
            using std::abs;
            assert(B1 != 0);
            assert(B2 != 0);
            using std::cerr;
            using std::endl;
            cerr << "AfterStep, i: " << i << endl;
            cerr << "AfterStep, j: " << j << endl;
            cerr << "AfterStep, B1 from q: " << q.get(i, j, 6) << endl;
            cerr << "AfterStep, B2 from q: " << q.get(i, j, 7) << endl;
            cerr << "AfterStep, B1: " << B1 << endl;
            cerr << "AfterStep, B2: " << B2 << endl;
            cerr << "AfterStep, A3 exact: " << A3_exact(angle, t, x, y) << endl;
            cerr << "AfterStep, A3: " << aux.get(i, j, 1) << endl;
//            assert(abs(B1) < 1.1);
//            assert(abs(B2) < 1.1);
//            assert(abs(B1 - q.get(i, j, 6)) < 0.2);
//            assert(abs(B2 - q.get(i, j, 7)) < 0.2);
            
            q.set(i, j, 6, B1);
            q.set(i, j, 7, B2);
        }
    } 
}
