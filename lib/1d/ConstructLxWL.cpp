#include <cmath>
#include "assert.h"            // for assert_eq.  Can be removed in future
#include "IniParams.h"
#include "IniParams.h"
#include "tensors.h"
#include "dog_math.h"

// Time-integrated right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = 0.
//
// This routine performs the Lax-Friedrich's flux
// splitting on a modified flux function, F.
void ConstructLxWL(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        const dTensorBC2& F,  // <-- new term: integrated flux
        dTensorBC2& Lstar,
        dTensorBC1& smax)
{

    // User supplied functions:
    void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void ProjectLeftEig( const dTensor1& Aux_ave, const dTensor1& Q_ave, 
        const dTensor2& Qvals, dTensor2& Wvals);
    void ProjectRightEig(const dTensor1& Aux_ave, const dTensor1& Q_ave, 
                         const dTensor2& Wvals, dTensor2& Qvals);
    void SetWaveSpd(const dTensor1& xedge, const dTensor1& Ql,
            const dTensor1& Qr, const dTensor1& Auxl, const dTensor1& Auxr,
            double& s1,double& s2);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);

    // Global wave speed
    void GlobalWaveSpd( const dTensorBC2& q, const dTensorBC2& aux, double& alpha1 );

    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Routine to deal with the silly mess where the Fluxes and the
    // Projections are all defined separately.
    void ConvertTranspose( const dTensor2& qin, dTensor2& qout );

    // Parameters for the current grid
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    // Determine size of WENO stencil
    const int ws = global_ini_params.get_space_order(); // Number of points for the weno-reconstruction
    const int r = (ws + 1) / 2;                 // order = 2*r-1
    assert_ge( mbc, r );

    // The flux, f_{i-1/2}.  Recall that the
    // flux lives at the nodal locations, i-1/2, so there is one more term in
    // that vector than on the original grid.
    dTensorBC2  fhat(mx+1, meqn, mbc );

    // Grid spacing
    const double   xlow = global_ini_params.get_xlow();
    const double     dx = global_ini_params.get_dx();

    double alpha1 = 0.;
    if( global_ini_params.get_global_alpha() )
    {
        // Global wave speed
        GlobalWaveSpd( q, aux, alpha1 );
    }

    // ---------------------------------------------------------
    // Compute fhat_{i-1/2}
    // ---------------------------------------------------------
#pragma omp parallel for
    for (int i= 1; i<= mx+1; i++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute "Roe" Averages.  We use simple arithmetic averages.
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,m) + q.get(i-1,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(iMax(maux, 1 ) );
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,ma) + aux.get(i-1,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( iMax(maux,1), ws+1         );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t(         ws+1, iMax(maux,1) );
        dTensor1 xvals( ws+1 );
        dTensor2 fvals( meqn, ws+1 );
        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int is = i-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( is )*dx - 0.5*dx;
            xvals.set( s, xi );
            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(is, m ) );
                fvals.set( m, s, F.get(is, m ) ); // <-- NEW part (sample integrated flux)
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(is, ma ) );
            }
        }

        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( Auxavg, Qavg, fvals, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------


        // -- Compute a local wave speed -- //

        dTensor1 xedge(1), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(iMax(1,maux)), Auxr(iMax(1,maux));
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i-1,m) );
            Qr.set(m, q.get(i  ,m) );
        }

        for( int m=1; m<= maux; m++)
        {
            Auxl.set(m, aux.get(i-1,m) );
            Auxr.set(m, aux.get(i  ,m) );
        }

        double s1,s2;
        SetWaveSpd(xedge, Ql, Qr, Auxl, Auxr, s1, s2);  // application specific
        const double alpha = Max( alpha1, Max( abs(s1), abs(s2) ) );
        smax.set( i, alpha  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here


        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s)      + l_alpha*wvals.get(m,s) )      );
            gm.set( m, s, 0.5*(gvals.get(m,ws-s+2) - l_alpha*wvals.get(m,ws-s+2) ) );
        }

        // --------------------------------------------------------------------
        // Part IV: Perform a WENO reconstruction on the characteristic vars.
        // --------------------------------------------------------------------
        dTensor2 dGp( meqn, 1 ), dGm( meqn, 1 );
        WenoReconstruct( gp, dGp );
        WenoReconstruct( gm, dGm );

        // add and convert back to the conserved quantities
        dTensor2 ghat( meqn, 1 );
        for (int m=1; m<=meqn; m++)
        {
            ghat.set(m, 1, dGp.get(m,1) + dGm.get(m,1) );
        }

        dTensor2 fhat_loc( ghat );
        ProjectRightEig(Auxavg, Qavg, ghat, fhat_loc);
        for( int m=1; m <= meqn; m++ )
        {
            fhat.set(i, m, fhat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_i = Lstar = -1/dx( fh_{i+1/2} - fh_{i-1/2} )
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {
        printf("Error: source-term not implemented for Lax-Wendroff method\n");
        exit(1);
        // Compute the source term.
//      SampleFunction( 1-mbc, mx+mbc, node, q, aux, Lstar, &SourceTermFunc);
//  #pragma omp parallel for
//          for (int i=1; i<=mx; i++)
//          {
//              for (int m=1; m<=meqn; m++)
//              {
//                  Lstar.set(i,m, Lstar.get(i,m) -(fhat.get(i+1,m)-fhat.get(i,m))/dx );
//              }
//          }
    }
    else
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        {
            for (int m=1; m<=meqn; m++)
            {
                Lstar.set(i,m, -(fhat.get(i+1,m)-fhat.get(i,m))/dx );
            }
        }
    }

    // ---------------------------------------------------------
    // Add extra contributions to Lstar
    // ---------------------------------------------------------
    // LstarExtra(aux,q,Lstar);

}
