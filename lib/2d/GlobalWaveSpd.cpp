#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "GlobalWaveSpd.h"

// Compute a global wave speed.
//
// alpha1 = max( |f'(q)| ), and
// alpha2 = max( |g'(q)| ), where the maximum is taken over each element q_{ij},
// and the hyperbolic problem is defined as
//
//     q_t + ( f(q) )_x + ( g(q) )_y = Psi.
//
// In order to speed up the calls to this function ( presumably by a factor of
// 2), one would need to redefine this function in a local application
// sub-directory.
//
// See also: SetWaveSpd for the local version of this function.
void GlobalWaveSpd(
    const dTensorBC3& q, 
    const dTensorBC3& aux, 
    double& alpha1, double& alpha2)
{

    // Grid and problem information
    const int mx     = dogParamsCart2.get_mx();
    const int my     = dogParamsCart2.get_my();
    const int meqn   = dogParams.get_meqn();
    const int maux   = dogParams.get_maux();
    const int mbc    = dogParamsCart2.get_mbc();

    // Needed to define derivatives
    const double dx    = dogParamsCart2.get_dx();
    const double dy    = dogParamsCart2.get_dy();
    const double xlow  = dogParamsCart2.get_xlow();
    const double ylow  = dogParamsCart2.get_ylow();

    alpha1 = alpha2 = 0.0;

    dTensor1 nvecx(2);
    nvecx.set(1, 1.0 ); nvecx.set(2, 0.0 );

    dTensor1 nvecy(2);
    nvecy.set(1, 0.0 ); nvecy.set(2, 1.0 );

    // TODO - replace this with a "paralell" version that uses locks each time
    // alpha1 and alpha2 are updated.  Right now, this is done in serial
// #pragma omp parallel for
    for( int i = 1; i <= mx; i++ )
    for( int j = 1; j <= my; j++ )
    {

        // Cell location
        dTensor1 xedge(2);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        // Solution values
        dTensor1 Ql(meqn), Qr(meqn);
        for( int me=1; me <= meqn; me++ )
        { Ql.set( me, q.get(i,j,me) ); }

        // Aux values
        dTensor1 Auxl(maux);
        for( int me=1; me <= maux; me++ )
        { Auxl.set( me, aux.get(i,j,me) ); }

        double s1 = 0.;
        double s2 = 0.;

        SetWaveSpd(nvecx, xedge, Ql, Ql, Auxl, Auxl, s1, s1 );
        alpha1 = Max( alpha1, s1 );

        SetWaveSpd(nvecy, xedge, Ql, Ql, Auxl, Auxl, s2, s2 );
        alpha2 = Max( alpha2, s2 );

    }

}
