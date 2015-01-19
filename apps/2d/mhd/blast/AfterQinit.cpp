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
    dTensorBC3& aux     = Qnew.ref_aux();
    //const double t      = Qnew.get_t();

    const int mx = global_ini_params.get_mx();
    const int my = global_ini_params.get_my();
    const double xlow = global_ini_params.get_xlow();
    //const double xhigh = global_ini_params.get_xhigh();
    const double ylow = global_ini_params.get_ylow();
    //const double yhigh = global_ini_params.get_yhigh();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

    for(int i = 1; i <= mx; ++i){
        for(int j = 1; j <= my; ++j){
            double x = xlow + (i - 0.5) * dx;
            double y = ylow + (j - 0.5) * dy;
            double tmp = 100.0 / sqrt(4.0*M_PI) * y;
            aux.set(i, j, 1,
                    tmp);

        }
    }
}
