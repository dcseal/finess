#include <cmath>
#include <cassert>
#include "tensors.h"
#include "StateVars.h"
#include "IniParams.h"


// This function is called ONCE per simulation, and immediately after 
// the initial conditions are set.
void AfterQinit( StateVars& Qnew )
{
    using std::exp;
    using std::pow;
    dTensorBC3& aux     = Qnew.ref_aux();
    
    const int mx = global_ini_params.get_mx();
    const int my = global_ini_params.get_my();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();
    const double mu = global_ini_params.get_mu();


    for(int i = 1; i <= mx; ++i){
        for(int j = 1; j <= my; ++j){
            double x = xlow + (i - 0.5) * dx;
            double y = ylow + (j - 0.5) * dy;
            aux.set(i, j, 1,
                    mu/(2.0*M_PI) * exp(0.5*(1.0-(x*x + y*y))) );
            
        }
    }

}
