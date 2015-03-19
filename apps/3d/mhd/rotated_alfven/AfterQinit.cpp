#include <cmath>
#include <cassert>
#include "tensors.h"
#include "StateVars.h"
#include "IniParams.h"


// This function is called ONCE per simulation, and immediately after 
// the initial conditions are set.
void AfterQinit( StateVars& Qnew )
{
    using std::cos;
    using std::sin;
    //dTensorBC3& q       = Qnew.ref_q();
    dTensorBC4& aux     = Qnew.ref_aux();
    //const double t      = Qnew.get_t();

    const int mx = global_ini_params.get_mx();
    const int my = global_ini_params.get_my();
    const int mz = global_ini_params.get_mz();
    const double xlow = global_ini_params.get_xlow();
    //const double xhigh = global_ini_params.get_xhigh();
    const double ylow = global_ini_params.get_ylow();
    //const double yhigh = global_ini_params.get_yhigh();
    const double zlow = global_ini_params.get_zlow();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();
    const double dz = global_ini_params.get_dz();
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};
    //    const double t[] = {-sin(phi), cos(phi), 0.0};
    //    const double r[] = {-cos(phi)*sin(theta), -sin(phi)*sin(theta), cos(theta)};


    const double b1av = n[0], b2av = n[1], b3av = n[2];
    for(int i = 1; i <= mx; ++i){
        for(int j = 1; j <= my; ++j){
            for(int k = 1; k <= mz; ++k){
                double x = xlow + (i - 0.5) * dx;
                double y = ylow + (j - 0.5) * dy;
                double z = zlow + (k - 0.5) * dz;
                double vnx = n[0] * x + n[1] * y + n[2] * z;
                aux.set(i, j, k, 1,
                        z * b2av - sin(phi)*sin(2.0*M_PI*vnx) / 20.0 / M_PI);
                aux.set(i, j, k, 2,
                        x * b3av + cos(phi)*sin(2.0*M_PI*vnx) / 20.0 / M_PI);
                aux.set(i, j, k, 3,
                        y * b1av + cos(2.0*M_PI*vnx) / 20.0 / M_PI / cos(theta) );

            }
        }
    }

}
