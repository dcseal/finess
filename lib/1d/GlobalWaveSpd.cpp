#include "dogdefs.h"
#include "dog_math.h"
#include "IniParams.h"
#include "IniParams.h"
#include "GlobalWaveSpd.h"

// Compute a global wave speed.
//
// alpha1 = max( |f'(q)| ), where the maximum is taken over each element q_{i},
// and the hyperbolic problem is defined as
//
//     q_t + ( f(q) )_x = Psi.
//
// In order to speed up the calls to this function ( presumably by a factor of
// 2), one would need to redefine this function in a local application
// sub-directory.
//
// See also: SetWaveSpd for the local version of this function.
void GlobalWaveSpd(
    const dTensorBC2& q, 
    const dTensorBC2& aux, 
    double& alpha1 )
{

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double xlow  = global_ini_params.get_xlow();

    alpha1 = 0.0;

    // TODO - replace this with a "paralell" version that uses locks each time
    // alpha1 and alpha2 are updated.  Right now, this is done in serial
// #pragma omp parallel for
    for( int i = 1; i <= mx; i++ )
    {

        // Cell location
        dTensor1 xedge(1);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );

        // Solution values
        dTensor1 Ql(meqn);
        // , Qr(meqn);
        for( int me=1; me <= meqn; me++ )
        { 
            Ql.set( me, q.get(i,me) ); 
//          Qr.set( me, q.get(i,me) ); 
        }

        // Aux values
        dTensor1 Auxl(maux), Auxr(maux);
        for( int me=1; me <= maux; me++ )
        { 
            Auxl.set( me, aux.get(i,me) ); 
//          Auxr.set( me, aux.get(i,me) ); 
        }

        double s1,s2;
        SetWaveSpd(xedge, Ql, Ql, Auxl, Auxl, s1, s2);  // application specific
        alpha1 = Max( alpha1, Max( fabs(s1), fabs(s2) ) );

    }

}
