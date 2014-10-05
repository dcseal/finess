#include <cmath>
#include <iostream>
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"

void AfterQinit( StateVars& Qnew );
{

    dTensorBC2& qnew    = Qnew.ref_q();
    dTensorBC2& aux     = Qnew.ref_aux();
    const double t      = Qnew.get_t();

    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();
    const int maux = aux.getsize(3);

    for(int i=1-mbc; i <= mx+mbc; i++ )
    for(int j=1-mbc; j <= my+mbc; j++ )
    {
        // flatten out any cells that are negative
        if( qnew.get(i,j,1,1) < 1e-13 )
        { 
            for( int k=1 ; k <= kmax; k++ )
            {
                qnew.set(i,j,1,k, 0.0 ); 
            }
        }
    }

//  void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);
//  if( global_ini_params.using_moment_limiter() )
//  { 
//      ApplyPosLimiter(aux, qnew); 
//  }
 

}
