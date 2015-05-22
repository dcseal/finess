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
    dTensorBC3& aux     = Qnew.ref_aux();
    //const double t      = Qnew.get_t();
    
    const int mx = global_ini_params.get_mx();
    const int my = global_ini_params.get_my();
    const int mbc = global_ini_params.get_mbc();
    const double xlow = global_ini_params.get_xlow();
    //const double xhigh = global_ini_params.get_xhigh();
    const double ylow = global_ini_params.get_ylow();
    //const double yhigh = global_ini_params.get_yhigh();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();
    const double angle = global_ini_params.get_angle();

    for(int i = 1-mbc; i <= mx+mbc; ++i){
        for(int j = 1-mbc; j <= my+mbc; ++j){
            double x = xlow + (i - 0.5) * dx;
            double y = ylow + (j - 0.5) * dy;
            double xi = x * cos(angle) + y * sin(angle);
            double eta = -x * sin(angle) + y * cos(angle);
            aux.set(i, j, 1,
                    xi <= 0 ?
                    0.75*eta - xi :
                    0.75*eta + xi );            
        }
    }

}
