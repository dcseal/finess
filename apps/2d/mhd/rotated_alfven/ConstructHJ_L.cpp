#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "ConstructHJ_L.h"

#include "StateVars.h"
#include "IniParams.h"
#include "assert.h"


#include <iostream>



static std::pair<double, double> max_speed_in_x_y_directions(const StateVars& Q);


// Reference: arxiv 1309.3344v2.  See equation (5.15).
//
void ConstructHJ_L(const StateVars& Q, dTensorBC3& Lauxstar)
{
    using std::pair;
    const dTensorBC3&   q = Q.const_ref_q();
    const dTensorBC3& aux = Q.const_ref_aux();


    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    assert(maux == 1);
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int    mbc = global_ini_params.get_mbc();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

    // Size of the WENO stencil
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );
    
    const pair<double, double> max_speed_x_y = max_speed_in_x_y_directions(Q);
    const double max_speed_x = max_speed_x_y.first;
    const double max_speed_y = max_speed_x_y.second;
    
    assert(maux == 1);
    assert(mbc >= 2 + r);
#pragma omp parallel for
    for(int i = 1 - 2 ; i <= mx + 2 ; ++i){
         dTensor2 weno_input(maux, ws);
         dTensor2 weno_result(maux, 1);
         for(int j = 1 - 2 ; j <= my + 2; ++j){
        
            // A^{3}_{,x}^{-}
            double A3xm;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i - r + k, j, 1) - aux.get(i - r + k - 1, j, 1)) / dx);
            WenoReconstruct(weno_input, weno_result);
            A3xm = weno_result.get(1, 1);

            // A^{3}_{,x}^{+}
            double A3xp;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i + r - (k - 1), j, 1) - aux.get(i + r - k, j, 1)) / dx);
            WenoReconstruct(weno_input, weno_result);
            A3xp = weno_result.get(1, 1);


            // A^{3}_{,y}^{-}
            double A3ym;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i, j - r + k, 1) - aux.get(i, j - r + k - 1, 1)) / dy);
            WenoReconstruct(weno_input, weno_result);
            A3ym = weno_result.get(1, 1);

            // A^{3}_{,y}^{+}
            double A3yp;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i, j + r - (k - 1), 1) - aux.get(i, j + r - k, 1)) / dy);
            WenoReconstruct(weno_input, weno_result);
            A3yp = weno_result.get(1, 1);
            
            double u1 = q.get(i, j, 2) / q.get(i, j, 1);
            double u2 = q.get(i, j, 3) / q.get(i, j, 1);
            Lauxstar.set(i, j, 1,
                    -u1 * 0.5 * (A3xm + A3xp) -u2 * 0.5 * (A3ym + A3yp)
                    +max_speed_x * 0.5 * (A3xp - A3xm)
                    +max_speed_x * 0.5 * (A3yp - A3ym)   );
        }
    }
// #pragma omp parallel for
//    for(int i = 1; i <= mx; ++i){
//        for(int j = 1; j <= my; ++j){
//            const double B1 = 
//                1.0 / (12.0 * dy)
//                * (aux.get(i, j - 2, 1) - 8*aux.get(i, j-1, 1) 
//                   + 8*aux.get(i, j+1, 1) - aux.get(i, j+2, 1));
//            const double B2 =
//                -1.0 / (12.0 * dx)
//                * (aux.get(i-2, j, 1) - 8*aux.get(i-1, j, 1)
//                   + 8*aux.get(i+1, j, 1) - aux.get(i+2, j, 1));
//
//            using std::abs;
//            assert(B1 != 0);
//            assert(B2 != 0);
//            using std::cerr;
//            using std::endl;
//            cerr << "ConstructHJ_L, B1: " << B1 << endl;
////            assert(abs(B1) < 1.1);
//            cerr << "ConstructHJ_L, B2: " << B2 << endl;
//            cerr << "ConstructHJ_L, B1 from q: " << q.get(i, j, 6) << endl;
//            cerr << "ConstructHJ_L, B2 from q: " << q.get(i, j, 7) << endl;
//            assert(abs(B2) < 2);
//        }
//    } 
    
    
}


static std::pair<double, double> max_speed_in_x_y_directions(const StateVars& Q){
    using std::abs;
    using std::max;
    using std::make_pair;

    const dTensorBC3& q = Q.const_ref_q();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    
    double max_x_speed = 0;
    double max_y_speed = 0;
    for(int i = 1; i <= mx; ++i){
        for(int j = 1; j <= my; ++j){
            assert(q.get(i, j, 1) != 0);
            double x_speed = std::abs(q.get(i, j, 2) / q.get(i, j, 1));
            double y_speed = std::abs(q.get(i, j, 3) / q.get(i, j, 1));
            max_x_speed = max(max_x_speed, x_speed);
            max_y_speed = max(max_y_speed, y_speed);
        }
    }
    using std::cerr;
    using std::endl;
    cerr << "max_x_speed: " << max_x_speed << endl;
    cerr << "max_y_speed: " << max_y_speed << endl;
    return make_pair(max_x_speed, max_y_speed);
}
