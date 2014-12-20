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
            double tmp = x < 0.05 ?
                -2.1826182 * x + 0.080921431 :
                -0.56418958 * x;
            aux.set(i, j, 1,
                    tmp);
            
        }
    }

}
