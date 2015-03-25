#include <cmath>
#include <cassert>
#include "tensors.h"
#include "StateVars.h"
#include "IniParams.h"


// This function is called ONCE per simulation, and immediately after 
// the initial conditions are set.
void AfterQinit( StateVars& Qnew )
{
    using std::sqrt;
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
    const int mbc = global_ini_params.get_mbc();

    for(int i = 1-mbc; i <= mx+mbc; ++i){
        for(int j = 1-mbc; j <= my+mbc; ++j){
            for(int k = 1-mbc; k <= mz+mbc; ++k){
                double x = xlow + (i - 0.5) * dx;
                double y = ylow + (j - 0.5) * dy;
                double z = zlow + (k - 0.5) * dz;
                double px, py, pz;

                px = 0.0;
                py = 0.0;
                pz = 100.0 / sqrt(4.0*M_PI) / sqrt(2) * y
                    -100.0 / sqrt(4.0*M_PI) / sqrt(2) * x;
                aux.set(i, j, k, 1, px);
                aux.set(i, j, k, 2, py);
                aux.set(i, j, k, 3, pz);

            }
        }
    }
}
