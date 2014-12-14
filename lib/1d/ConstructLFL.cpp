#include <cmath>
#include "assert.h"            // for assert_eq.  Can be removed in future
#include "IniParams.h"
#include "IniParams.h"
#include "tensors.h"
#include "dog_math.h"
#include "ConstructL.h"
#include "StateVars.h"

// Construct flux values for a Lax-Friedrich's flux
//
// fhat(1:mx+1, meqn, mbc), fhat(i) = fhat_{i-1/2}.
//
void ConstructLFL( const StateVars& Q, dTensorBC2& fhat )
{

    const dTensorBC2&    q = Q.const_ref_q  ();
    const dTensorBC2&  aux = Q.const_ref_aux();

    // Parameters for the current grid (could also use dogParams here)
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    // Grid spacing
    const double   xlow = global_ini_params.get_xlow();
    const double     dx = global_ini_params.get_dx();

    // (Global) Wave speed
    double alpha = 0.;
    GlobalWaveSpd( q, aux, alpha );

    // ---------------------------------------------------------
    // Compute fhat_{i-1/2}
    // ---------------------------------------------------------
#pragma omp parallel for
    for (int i= 1; i<= mx+1; i++)
    {

        // Two grid points:
        dTensor1 xvals( 2 );
        xvals.set(1, xlow+double(i)*dx - 0.5*dx );
        xvals.set(2, xlow+double(i)*dx + 0.5*dx );

        // Sample q and aux over the stencil:
        dTensor2 qvals_t( 2, meqn ), auxvals_t( 2, maux );
        for( int m=1; m <= meqn; m++ )
        {
            qvals_t.set( 1, m, q.get(i-1, m ) );
            qvals_t.set( 2, m, q.get(i  , m ) );
        }
        for( int ma=1; ma <= maux; ma++ )
        {
            auxvals_t.set( 1, ma, aux.get(i-1, ma ) );
            auxvals_t.set( 2, ma, aux.get(i  , ma ) );
        }

        // Sample f over the stencil:
        dTensor2 fvals_t( 2, meqn );
        FluxFunc( xvals, qvals_t, auxvals_t, fvals_t );

        // Save the flux function
        for( int m=1; m <= meqn; m++ )
        {
            fhat.set( i, m, 0.5*( (fvals_t.get(1,m)+fvals_t.get(2,m)) - alpha*(q.get(i,m) - q.get(i-1,m))) );
        }
    }

}
