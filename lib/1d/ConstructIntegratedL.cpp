#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "assert.h"

void Diff1( double dx, const dTensor2& f, dTensor1& fx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils:
    assert_eq( f.getsize( 2 ), 6 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = (  f.get( m, 1 ) - f.get( m, 6 ) )*(1.0/12.0);
        tmp       += ( -f.get( m, 2 ) + f.get( m, 4 ) )*(2.0/ 3.0);
        fx.set( m, tmp / dx );
    }

}

void Diff2( double dx, const dTensor2& f, dTensor1& fxx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils:
    assert_eq( f.getsize( 2 ), 6 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = ( -f.get( m, 1 ) - f.get( m, 6 ) )*(1.0/12.0);
        tmp       += (  f.get( m, 2 ) + f.get( m, 4 ) )*(4.0/ 3.0);
        tmp       += ( -f.get( m, 3 )                 )*(5.0/ 2.0);
        fxx.set( m, tmp / (dx*dx) );
    }

}

void ConstructIntegratedL( double dt, const dTensor2& node, 
    const dTensorBC2& aux, const dTensorBC2& q,                       // current state
    dTensorBC2& f, dTensorBC2& fx, dTensorBC2& fxx, dTensorBC2& qx,   // recycled storage
    dTensorBC1& smax, dTensorBC2& F)
{


    void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void SampleFunction( 
        int istart, int iend,
        const dTensor2& node,
        const dTensorBC2& qin, 
        const dTensorBC2& auxin,  
        dTensorBC2& Fout,
        void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

    // Problem dimensions (TODO - the boundary data either a) needs one more
    // point, or b) needs to double the number of ghost cells )
    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int mbc    = q.getmbc();

    const double dx  = dogParamsCart1.get_dx();

    // Sample the flux function on the entire domain:
    SampleFunction( 1-mbc, mx+mbc, node, q, aux, f, &FluxFunc );

// TODO  - allow for different sized stencils
const int      mpts_sten = 2*mbc;  assert_eq( mpts_sten, 6 );
const int half_mpts_sten =   mbc;  assert_eq( mpts_sten, 3 );

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1; i <= mx; i++ )
    {
        // Save the flux function:
        dTensor2 Fvals( meqn, mpts_sten );
        dTensor2 qvals( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                Fvals.set( m, r, f.get( i+s, m ) );
                qvals.set( m, r, q.get( i+s, m ) );
                s++;
            }
        }

        // Storage for local derivatives:
        dTensor1 fx_val  ( meqn );
        dTensor1 qx_val  ( meqn );
        dTensor1 fxx_val ( meqn );

        // Compute a FD approximation to the derivatives:
        Diff1( dx, Fvals, fx_val  );
        Diff1( dx, qvals, qx_val  );
        Diff2( dx, Fvals, fxx_val );

        // Save the derivatives for later use:
        for( int m=1; m <= meqn; m++ )
        {
            fx.set  ( i, m,  fx_val.get(m) );
            qx.set  ( i, m,  qx_val.get(m) );
            fxx.set ( i, m, fxx_val.get(m) );
        }

    }

    // TODO - something needs to be done about the boundary data!!!
printf("Error: ConstructIntegratedL has not been finished\n");
exit(1);

}
