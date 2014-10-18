#include <cmath>
#include <cassert>
#include <iostream>
#include "tensors.h"
#include "StateVars.h"
#include "IniParams.h"

#include "app_util.h"


// Function that is called before each stage
void BeforeStep(double dt, StateVars& Q)
{
//    using std::cos;
//    using std::sin;
//    using std::abs;
//    using std::cerr;
//    using std::endl;
//
//    const dTensorBC3& q       = Q.const_ref_q();
//    dTensorBC3& aux     = Q.ref_aux();
//    //const double t      = Qnew.get_t();
//    
//    double t = Q.get_t();
//    cerr << "Stage at time:  " << t << endl;
//
//
//    const int mx = global_ini_params.get_mx();
//    const int my = global_ini_params.get_my();
//    const double xlow = global_ini_params.get_xlow();
//    //const double xhigh = global_ini_params.get_xhigh();
//    const double ylow = global_ini_params.get_ylow();
//    //const double yhigh = global_ini_params.get_yhigh();
//    const double dx = global_ini_params.get_dx();
//    const double dy = global_ini_params.get_dy();
//    const double angle = global_ini_params.get_angle();
//    const int check_padding = 0; 
    

//    cerr << "BeforeStep,  A[0, 0]: " << aux.get(0, 0, 1) << endl;
//    cerr << "BeforeStep,  A[0, 1]: " << aux.get(0, 0, 1) << endl;

//#pragma omp parallel for
//    for(int i = 1 + check_padding; i <= mx - check_padding; ++i){
//        for(int j = 1 + check_padding; j <= my - check_padding; ++j){
//            const double B1 = 
//                1.0 / (12.0 * dy)
//                * (aux.get(i, j - 2, 1) - 8*aux.get(i, j-1, 1) 
//                   + 8*aux.get(i, j+1, 1) - aux.get(i, j+2, 1));
//            const double B2 =
//                -1.0 / (12.0 * dx)
//                * (aux.get(i-2, j, 1) - 8*aux.get(i-1, j, 1)
//                   + 8*aux.get(i+1, j, 1) - aux.get(i+2, j, 1));
//
//            assert(B1 != 0);
//            assert(B2 != 0);
//            cerr << "BeforeStep, i: " << i << endl;
//            cerr << "BeforeStep, j: " << j << endl;
//            cerr << "BeforeStep A3: " << aux.get(i, j, 1) << endl;
//            double x, y;
//            ij_to_xy(xlow, ylow, dx, dy, i, j, x, y);
//            double a3e = A3_exact(angle, t, x, y);
//            cerr << "BeforeStep A3 exact: "
//                 << a3e
//                 << endl;
//            assert(abs(aux.get(i, j, 1) - a3e) < 0.1);
//
//            cerr << "BeforeStep, B1: " << B1 << endl;
////            assert(abs(B1) < 1.1);
//            cerr << "BeforeStep, B2: " << B2 << endl;
//
//            cerr << "BeforeStep, B1 from q: " << q.get(i, j, 6) << endl;
//            cerr << "BeforeStep, B2 from q: " << q.get(i, j, 7) << endl;
//            assert(abs(B1 - q.get(i, j, 6))< 0.1);
//            assert(abs(B2 - q.get(i, j, 7)) < 0.1);
//
//        }
//    } 
 }
