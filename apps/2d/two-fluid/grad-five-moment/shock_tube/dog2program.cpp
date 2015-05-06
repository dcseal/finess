#include "dogdefs.h"
#include <string.h>          // For strcpy and strcat (some compilers don't need this)
#include "IniParams.h"
#include "StateVars.h"


// Function that is called after initial condition
void AfterQinit(StateVars& Qstate)
{

  //DO NOTHING


}
#include <fstream>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor2& xpts, 
	     dTensor2& auxvals)
{
  //DO NOTHING
}
#include "Components.h"         // Easier index into components of solution
#include "tensors.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     2D five-moment two-fluid equations
//
// The flux can be decoupled into a total of "three" smaller fluxes:
//
//      a) the flux for the electrons
//      b) the flux for the ions
//      c) the flux for Maxwell's equations
//
void FluxFunc(const dTensor2& xpts,
	      const dTensor2& Q,
	      const dTensor2& Aux,
	      dTensor3& flux)
{

    void FiveMomentFluxFunc( int n_offset, const dTensor2& Q, dTensor3& flux);
    FiveMomentFluxFunc(0, Q, flux);
    FiveMomentFluxFunc(5, Q, flux);

    void MaxwellFluxFunc( int n_offset, const dTensor2& Q, dTensor3& flux);
    MaxwellFluxFunc(10, Q, flux);

    // TODO - I'm not sure what this part is about ... -DS
    void AdvectionFluxFunc( const dTensor2& Q, dTensor3& flux, int advIdx, int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
    AdvectionFluxFunc(Q, flux, _entropy_i, _rho_i);
    AdvectionFluxFunc(Q, flux, _entropy_e, _rho_e);

}

void FluxFunc1( const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor2& flux)
{

    void FiveMomentFluxFunc1( int n_offset, const dTensor2& Q, dTensor2& flux);
    FiveMomentFluxFunc1(0, Q, flux);
    FiveMomentFluxFunc1(5, Q, flux);

    void MaxwellFluxFunc1( int n_offset, const dTensor2& Q, dTensor2& flux);
    MaxwellFluxFunc1(10, Q, flux);

    void AdvectionFluxFunc1( const dTensor2& Q, dTensor2& flux, int advIdx, int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
    AdvectionFluxFunc1(Q, flux, _entropy_i, _rho_i);
    if(Q.getsize(2)<_entropy_e) return;
    AdvectionFluxFunc1(Q, flux, _entropy_e, _rho_e);

}

void FluxFunc2( const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor2& flux)
{

    void FiveMomentFluxFunc2( int n_offset, const dTensor2& Q, dTensor2& flux);
    FiveMomentFluxFunc2(0, Q, flux);
    FiveMomentFluxFunc2(5, Q, flux);

    void MaxwellFluxFunc2( int n_offset, const dTensor2& Q, dTensor2& flux);
    MaxwellFluxFunc2(10, Q, flux);

    void AdvectionFluxFunc2( const dTensor2& Q, dTensor2& flux, int advIdx, int rhoIdx);
    if(Q.getsize(2)<_entropy_i) return;
    AdvectionFluxFunc2(Q, flux, _entropy_i, _rho_i);
    if(Q.getsize(2)<_entropy_e) return;
    AdvectionFluxFunc2(Q, flux, _entropy_e, _rho_e);

}

#include "Components.h"
// This is a user-supplied routine that projects
// Q onto the left eigenvectors of the flux 
// Jacobian; the result is stored in W
//
class dTensor1;
class dTensor2;
void ProjectLeftEig( int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{    
    void ProjectLeftEig_FiveMoment( int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
    ProjectLeftEig_FiveMoment(ixy, 0, Q_ave, Q, W);
    ProjectLeftEig_FiveMoment(ixy, 5, Q_ave, Q, W);

    void ProjectLeftEig_Maxwell(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
    ProjectLeftEig_Maxwell(ixy, 10, Q_ave, Q, W);

    void ProjectLeftConvectedScalar(int idx,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
    ProjectLeftConvectedScalar(_entropy_i, Q_ave, Q, W);
    ProjectLeftConvectedScalar(_entropy_e, Q_ave, Q, W);
}
#include "Components.h"
// This is a user-supplied routine that projects
// W onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Q
//
class dTensor1;
class dTensor2;
void ProjectRightEig( int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{    
    void ProjectRightEig_FiveMoment(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
    ProjectRightEig_FiveMoment(ixy, 0, Q_ave, W, Q);
    ProjectRightEig_FiveMoment(ixy, 5, Q_ave, W, Q);

    void ProjectRightEig_Maxwell(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
    ProjectRightEig_Maxwell(ixy, 10, Q_ave, W, Q);

    void ProjectRightConvectedScalar(int idx,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
    ProjectRightConvectedScalar(_entropy_i, Q_ave, W, Q);
    ProjectRightConvectedScalar(_entropy_e, Q_ave, W, Q);
}
#include <cmath>
#include "dog_math.h"
#include "tensors.h"

#include "IniParams.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
// Ideal MHD equations
void SetWaveSpd(const dTensor1& nvec, 
		const dTensor1& xedge,
		const dTensor1& Ql, 
		const dTensor1& Qr,
		const dTensor1& Auxl, 
		const dTensor1& Auxr,
		double& s1,double& s2)
{

  // Gas constant
  const double gamma = global_ini_params.get_gamma();
    
  // Left states
  const double rhol = Ql.get(1);
  const double u1l  = Ql.get(2)/rhol;
  const double u2l  = Ql.get(3)/rhol;
  const double u3l  = Ql.get(4)/rhol;
  const double energyl = Ql.get(5);
  const double B1l  = Ql.get(6);
  const double B2l  = Ql.get(7);
  const double B3l  = Ql.get(8);
  const double Bm2l = B1l*B1l + B2l*B2l + B3l*B3l;
  const double um2l = u1l*u1l + u2l*u2l + u3l*u3l;
  const double pressl = (gamma-1.0e0)*(energyl  -0.5*rhol*um2l - 0.5*Bm2l);
  
  // Right states
  const double rhor = Qr.get(1);
  const double u1r  = Qr.get(2)/rhor;
  const double u2r  = Qr.get(3)/rhor;
  const double u3r  = Qr.get(4)/rhor;
  const double energyr = Qr.get(5);
  const double B1r  = Qr.get(6);
  const double B2r  = Qr.get(7);
  const double B3r  = Qr.get(8);
  const double Bm2r = B1r*B1r + B2r*B2r + B3r*B3r;
  const double um2r = u1r*u1r + u2r*u2r + u3r*u3r;
  const double pressr = (gamma-1.0e0)*(energyr - 0.5*rhor*um2r - 0.5*Bm2r);
  
  // Average states
  const double rho    = 0.5e0*(rhol+rhor);
  const double u1     = 0.5e0*(u1l+u1r);
  const double u2     = 0.5e0*(u2l+u2r);
  const double u3     = 0.5e0*(u3l+u3r);
  const double press  = 0.5e0*(pressl+pressr);
  const double B1     = 0.5e0*(B1l+B1r);
  const double B2     = 0.5e0*(B2l+B2r);
  const double B3     = 0.5e0*(B3l+B3r);
  const double Bm2    = B1*B1 + B2*B2 + B3*B3;
  
  // Sound speed squared
  const double a2l = fabs(gamma*pressl/rhol);
  const double a2r = fabs(gamma*pressr/rhor);
  const double a2  = fabs(gamma*press/rho);

  // Sound speeds
  const double nmag = sqrt(pow(nvec.get(1),2) + pow(nvec.get(2),2));
  //  const double cl = nmag*sqrt(fabs(gamma*pressl/rhol));
  //  const double cr = nmag*sqrt(fabs(gamma*pressr/rhor));
  //  const double c  = nmag*sqrt(fabs(gamma*press/rho));

  // Fast magnetosonic speeds
  const double cf_l = nmag*sqrt(0.5*(a2l + Bm2l/rhol + 
		   sqrt(fabs(pow(a2l+Bm2l/rhol,2) - 4.0*a2l*B1l*B1l/rhol))));
  const double cf_r = nmag*sqrt(0.5*(a2r + Bm2r/rhor + 
		   sqrt(fabs(pow(a2r+Bm2r/rhor,2) - 4.0*a2r*B1r*B1r/rhor))));
  const double cf   = nmag*sqrt(0.5*(a2 + Bm2/rho + 
		   sqrt(fabs(pow(a2+Bm2/rho,2) - 4.0*a2*B1*B1/rho))));
  

  
  // normal velocities
  const double un  = nvec.get(1)*u1  + nvec.get(2)*u2;
  const double unl = nvec.get(1)*u1l + nvec.get(2)*u2l;
  const double unr = nvec.get(1)*u1r + nvec.get(2)*u2r;

  // Minimum speed
  s1 = Min(unl-cf_l, un-cf);
  
  // Maximum speed
  s2 = Max(unr+cf_r, un+cf);
}
#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// Function that is called after a full time step (i.e., after all stages are complete)
void AfterFullTimeStep(double dt, StateVars& Qold, StateVars& Q)
{
//  const int   mx   = q.getsize(1);
//  const int   my   = q.getsize(2);
//  const int meqn   = q.getsize(3);
//  const int mbc    = q.getmbc();
//  const int maux   = aux.getsize(3);
}
#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// Function that is called after each stage
void AfterStep(double dt, StateVars& Q )
{

//  const int     mx = q.getsize(1);
//  const int     my = q.getsize(2);
//  const int   meqn = q.getsize(3);
//  const int   maux = aux.getsize(3);  
}
#include <cmath>
#include<iostream>
#include "stdio.h"
#include "dog_math.h"
#include "IniParams.h"
#include "constants.h"
#include "tensors.h"
using namespace std;

// MPP limiter.  This needs to be written individually for each application
// that attempts to use it.
//
// This limiter works by considering a high order flux, (fHat, gHat), and
// low-order flux, (fLF, gLF), and then creates a new flux of the form:
//
//     fhat = theta fhat + (1-theta) fLF
//     ghat = theta ghat + (1-theta) gLg
//
// so that the conservative update, 
//
// q = q - dt/dx( f_{i+1/2} - f_{i-1/2} ) - dt/dy( g_{j+1/2} - g_{j-1/2} )
//
// is conservative.
//
// See: "http://arxiv.org/abs/1411.0328" and references therein for more
// details.
void ApplyMPPLimiter2D( 
        const double dt, const dTensorBC3& q, 
        const dTensorBC3& fLF, const dTensorBC3& gLF,
        dTensorBC3& fHat, dTensorBC3& gHat )
{
    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();

    printf("WARNING: no positivity limiter has been implemented for this application\n");

}
#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
//
// Function that is called before a full time step
void BeforeFullTimeStep(double dt, 
		       StateVars& Qold, 
		       StateVars& Q)
{

//  const int   mx   = q.getsize(1);
//  const int   my   = q.getsize(2);
//  const int meqn   = q.getsize(3);
//  const int mbc    = q.getmbc();
//  const int maux   = aux.getsize(3);

}
#include "tensors.h"
#include "StateVars.h"

// *TEMPLATE*
// 
// Function that is called before each stage
void BeforeStep(double dt, StateVars& Q)
{
//  const int     mx = q.getsize(1);
//  const int     my = q.getsize(2);
//  const int   meqn = q.getsize(3);
//  const int   maux = aux.getsize(3);  
}
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "StateVars.h"

#include "IniParams.h"
using namespace std;

void ConSoln( const StateVars& Q )
{

    const dTensorBC3&    q = Q.const_ref_q  ();
    const dTensorBC3&  aux = Q.const_ref_aux();
    const double         t = Q.get_t();

    string outputdir = global_ini_params.get_output_dir();

    // Size of the solution
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();

    // Grid information:
    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();

    string fname1 = outputdir+"/conservation.dat";
    ofstream write_file1,write_file2;
    dTensor1 qsum(meqn);
    dTensor1 res_sum(meqn);

    if( t==0 ) 
    {
        write_file1.open(fname1.c_str(), ofstream::out);
    }
    else
    {
        write_file1.open(fname1.c_str(), ofstream::app);
    }

    // -----------------
    // CONSERVATION
    // -----------------
    if( global_ini_params.get_mcapa() < 1 ) // without capacity function
    {
        for (int m=1; m<=meqn; m++)
        {

            qsum.set(m,0.0);
            for (int i=1; i<=mx; i++)
            for (int j=1; j<=my; j++)
            {
                const double x    = xlow + (double(i)-0.5)*dx;
                const double y    = ylow + (double(j)-0.5)*dy;
                const double qtmp = q.get(i,j,m);
                qsum.set(m, qsum.get(m) + dx*dy*qtmp );
            }
        }
    }
    else // with capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m, 0.0);

            for (int i=1; i<=mx; i++)
            for (int j=1; j<=my; j++)
            {
                const double x  = xlow + (double(i)-0.5)*dx;
                const double y  = ylow + (double(j)-0.5)*dy;

                double qtmp = q.get(i,j,m);
                double atmp = aux.get(i,j, global_ini_params.get_mcapa() );
                qsum.set(m, (qsum.get(m) + atmp*dx*dy*qtmp) );
            }
        }
    }

    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        if (abs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
        write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
    write_file1 << endl;

}
#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "assert.h"
#include "IniParams.h"
#include "StateVars.h"

// Central difference formulae
// see $FINESS/lib/WenoReconstruct.cpp
void Diff1( double dx, const dTensor2& f, dTensor1& fx );
void Diff2( double dx, const dTensor2& f, dTensor1& fxx );
double Diff1( double dx, double f1, double f2, double f3, double f4, double f5 );

// User supplied functions defining the Flux function, Jacobian, and
// Hessian of the flux function.
void FluxFunc(const dTensor2&,const dTensor2&,const dTensor2&,dTensor3&);
void DFluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux,
    dTensor4& Dflux);
void D2FluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux,
    dTensor5& D2flux);

// Used for construcing the flux function
void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC4& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&));


// This function computes the (linear) finite difference approximation to the
// integrated flux function on the conserved variables.  It requires knowledge
// of FluxFunc (1st-order), DFluxFunc (2nd-order) and D2FluxFunc (3rd-order)
// in order to compute the expansion:
//
//     F := f - dt/2 ( A (f_x + g_y ) )
//            + \cdots.
//
//     G := g - dt/2 ( B (f_x + g_y ) )
//            + \cdots.
//
// Where the flux Jacobian and Hessian are defined as:
//
//      A := \partial   f / \partial q,   and 
//      B := \partial   g / \partial q,   and 
//    A_q := \partial^2 f / \partial^2 q, and
//    B_q := \partial^2 g / \partial^2 q.
//
// Higher order methods would require further terms to be defined here,
// including "super"-Hessians.
//
// Lax-Friedrichs flux splitting + WENO reconstruction can then be applied 
// to the integrated flux function, F and G to define an update of the form:
//
//     q^{n+1} = q^n - \dt ( F_x + G_y )
//
// See also: DFluxFunc and D2FluxFunc.
void ConstructIntegratedR( double dt, const StateVars& Q,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G)
{

    const dTensorBC3& q   = Q.const_ref_q  ();
    const dTensorBC3& aux = Q.const_ref_aux();

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double dy    = global_ini_params.get_dy();
    const double xlow  = global_ini_params.get_xlow();
    const double ylow  = global_ini_params.get_ylow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC4 R( mx, my, meqn, 2, mbc );  // place-holder for the flux function
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, R, &FluxFunc );

// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );

const int ndim = 2;

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    for( int j = 1-mbc_small; j <= my+mbc_small; j++ )
    {

        // Physical location for this current value:
        dTensor2 xpts( 1, ndim );
        xpts.set( 1, 1, xlow + double(i)*dx - 0.5*dx );
        xpts.set( 1, 2, ylow + double(j)*dy - 0.5*dy );

        // Save the flux function:
        dTensor2 Fvals  ( meqn, mpts_sten );
        dTensor2 Gvals  ( meqn, mpts_sten );
        dTensor2 qvalsx ( meqn, mpts_sten );
        dTensor2 qvalsy ( meqn, mpts_sten );

        for( int m=1; m <= meqn; m++ )
        {
            int s = -half_mpts_sten+1;
            for( int r = 1; r <= mpts_sten; r++ )
            {
                Fvals.set( m, r, R.get( i+s, j, m, 1 ) ); // for computing f_x
                Gvals.set( m, r, R.get( i, j+s, m, 2 ) ); // for computing g_y
                qvalsx.set( m, r, q.get( i+s, j, m ) );
                qvalsy.set( m, r, q.get( i, j+s, m ) );
                s++;
            }
        }

        // Storage for local derivatives:
        dTensor1 fx_val  ( meqn ), qx_val( meqn );
        dTensor1 gy_val  ( meqn ), qy_val( meqn );

        // Compute a FD approximation to the derivatives:
        Diff1( dx, Fvals,  fx_val  );
        Diff1( dx, qvalsx, qx_val  );

        Diff1( dy, Gvals,  gy_val  );
        Diff1( dy, qvalsy, qy_val  );

        // Construct the first product: A f_x:
        dTensor4 A( 1, meqn, meqn, 2 );             // ( f'(q), g'(q) )
        dTensor2 q_transpose( 1, meqn );
        dTensor2 a_transpose( 1, maux );

        for( int m=1; m <= meqn; m++ )
        {
            q_transpose.set(1, m, q.get(i,j,m) );
        }
        for( int m=1; m <= maux; m++ )
        {
            a_transpose.set(1, m, aux.get(i,j,m) );
        }

        // Compute the Jacobian:
        DFluxFunc(xpts, q_transpose, a_transpose, A);

        // Compute the product: f'(q)*(f_x+g_y) + g'(q)*(f_x+g_y)
        dTensor1 f_t( meqn ), g_t( meqn );
        for( int m1=1; m1 <= meqn; m1++ )
        {
            double tmp1 = 0.;
            double tmp2 = 0.;
            for( int m2=1; m2 <= meqn; m2++ )
            {
                tmp1 += -(A.get(1, m1, m2, 1)) * ( fx_val.get(m2) + gy_val.get(m2));
                tmp2 += -(A.get(1, m1, m2, 2)) * ( fx_val.get(m2) + gy_val.get(m2));
            }
            f_t.set( m1, tmp1 );
            g_t.set( m1, tmp2 );
        }

        // ---  Third-order terms --- //
        dTensor1 f_tt( meqn );   f_tt.setall(0.);
        dTensor1 g_tt( meqn );   g_tt.setall(0.);
        if( global_ini_params.get_time_order() > 2 )
        {

            // Hessian
            dTensor5 H( 1, meqn, meqn, meqn, 2 );
            D2FluxFunc(xpts, q_transpose, a_transpose, H);

            // ----------------------------------- //
            // Part I: Compute (f_x + g_y)_{,t}
            // ----------------------------------- //

            dTensor1 fx_plus_gy_t( meqn ); fx_plus_gy_t.setall(0.);

            // Start with a term that get's used frequently:
            dTensor1 fx_plus_gy( meqn );
            for( int m =1; m <= meqn; m++ )
            {
                double tmp = fx_val.get(m) + gy_val.get(m);
                fx_plus_gy.set(m, tmp );
            }

            // Second-derivatives:
            dTensor1 fxx_val( meqn ), gyy_val( meqn );
            Diff2( dx, Fvals,  fxx_val  );
            Diff2( dy, Gvals,  gyy_val  );

            // Cross - derivaties
            dTensor1 fxy_val  ( meqn );  fxy_val.setall(0.);
            dTensor1 gxy_val  ( meqn );  gxy_val.setall(0.);

            // 2nd-order stencil (for mixed derivatives)
            for( int m=1; m <= meqn; m++ )
            {
                double tmp1  = 0.5*(R.get(i+1,j+1,m,1)-R.get(i-1,j+1,m,1))/dx;
                       tmp1 -= 0.5*(R.get(i+1,j-1,m,1)-R.get(i-1,j-1,m,1))/dx;
                       tmp1 *= 0.5/dy;
                fxy_val.set(m, tmp1);

                double tmp2  = 0.5*(R.get(i+1,j+1,m,2)-R.get(i-1,j+1,m,2))/dx;
                       tmp2 -= 0.5*(R.get(i+1,j-1,m,2)-R.get(i-1,j-1,m,2))/dx;
                       tmp2 *= 0.5/dy;
                gxy_val.set(m, tmp2);
            }

            // Compute terms that get multiplied by 
            //     \pd2{ f }{ q } and \pd2{ g }{ q }.
            for( int m =1; m <= meqn; m++ )
            {
                double tmp = 0.;

                // Terms that get multiplied by the Hessian:
                for( int m1=1; m1 <= meqn; m1++ )
                for( int m2=1; m2 <= meqn; m2++ )
                {

                    tmp += H.get(1,m,m1,m2,1)*qx_val.get(m1)*fx_plus_gy.get(m2);
                    tmp += H.get(1,m,m1,m2,2)*qy_val.get(m1)*fx_plus_gy.get(m2);
                }

                // Terms that get multiplied by f'(q) and g'(q):
                for( int m1=1; m1 <= meqn; m1++ )
                {

                    tmp += A.get(1,m,m1,1)*( fxx_val.get(m1)+gxy_val.get(m1) );
                    tmp += A.get(1,m,m1,2)*( fxy_val.get(m1)+gyy_val.get(m1) );
                }

                fx_plus_gy_t.set( m, tmp );
            }


            // ----------------------------------- //
            // Part II: Compute 
            //      f'(q) * fx_plus_gy_t and 
            //      g'(q) * fx_plus_gy_t
            // ----------------------------------- //

            // Add in the third term that gets multiplied by A:
            for( int m1=1; m1 <= meqn; m1++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += A.get(1,m1,m2,1)*fx_plus_gy_t.get(m2);
                    tmp2 += A.get(1,m1,m2,2)*fx_plus_gy_t.get(m2);
                }
                f_tt.set( m1, tmp1 );
                g_tt.set( m1, tmp2 );
            }

            // ----------------------------------------------- //
            // Part III: Add in contributions from
            //      f''(q) * (fx_plus_gy, fx_plus_gy ) and 
            //      g''(q) * (fx_plus_gy, fx_plus_gy ).
            // ----------------------------------------------- //
            for( int m =1; m <= meqn; m++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;

                // Terms that get multiplied by the Hessian:
                for( int m1=1; m1 <= meqn; m1++ )
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += H.get(1,m,m1,m2,1)*fx_plus_gy.get(m1)*fx_plus_gy.get(m2);
                    tmp2 += H.get(1,m,m1,m2,2)*fx_plus_gy.get(m1)*fx_plus_gy.get(m2);
                }

                f_tt.set( m, f_tt.get(m) + tmp1 );
                g_tt.set( m, g_tt.get(m) + tmp2 );
            }

        }

        // FINAL STEP: save the time-integrated values:
        for( int m=1; m<=meqn; m++ )
        {
            F.set(i,j,m, R.get(i,j,m,1) + 0.5*dt*(f_t.get(m) + dt/3.0*f_tt.get(m)) );
            G.set(i,j,m, R.get(i,j,m,2) + 0.5*dt*(g_t.get(m) + dt/3.0*g_tt.get(m)) );
        }


    }

}

void LocalIntegrate( 
    int nterms, double dx, double dy, double xc, double yc,
    int meqn, int maux, int mpts_sten, int half_mpts_sten,
    const int i, const int j, const dTensorBC3& q, const dTensorBC3& aux, 
    const dTensorBC4& R, 
    dTensor1& f_t, dTensor1& f_tt,
    dTensor1& g_t, dTensor1& g_tt
    )
{

    // Problem dimension (used for setting xpts)
    const int ndim = 2;

    // Physical location for this current value:
    dTensor2 xpts( 1, ndim );
    xpts.set( 1, 1, xc );
    xpts.set( 1, 2, yc );

    // Save the flux function:
    dTensor2 Fvals  ( meqn, mpts_sten );
    dTensor2 Gvals  ( meqn, mpts_sten );
    dTensor2 qvalsx ( meqn, mpts_sten );
    dTensor2 qvalsy ( meqn, mpts_sten );

    for( int m=1; m <= meqn; m++ )
    {
        int s = -half_mpts_sten+1;
        for( int r = 1; r <= mpts_sten; r++ )
        {
            Fvals.set( m, r, R.get( i+s, j, m, 1 ) ); // for computing f_x
            Gvals.set( m, r, R.get( i, j+s, m, 2 ) ); // for computing g_y
            qvalsx.set( m, r, q.get( i+s, j, m ) );
            qvalsy.set( m, r, q.get( i, j+s, m ) );
            s++;
        }
    }

    // Storage for local derivatives:
    dTensor1 fx_val  ( meqn ), qx_val( meqn );
    dTensor1 gy_val  ( meqn ), qy_val( meqn );

    // Compute a FD approximation to the derivatives:
    Diff1( dx, Fvals,  fx_val  );
    Diff1( dx, qvalsx, qx_val  );

    Diff1( dy, Gvals,  gy_val  );
    Diff1( dy, qvalsy, qy_val  );

    // Construct the first product: A f_x:
    dTensor4 A( 1, meqn, meqn, 2 );             // ( f'(q), g'(q) )
    dTensor2 q_transpose( 1, meqn );
    dTensor2 a_transpose( 1, maux );

    for( int m=1; m <= meqn; m++ )
    {
        q_transpose.set(1, m, q.get(i,j,m) );
    }
    for( int m=1; m <= maux; m++ )
    {
        a_transpose.set(1, m, aux.get(i,j,m) );
    }

    // Compute the Jacobian:
    DFluxFunc(xpts, q_transpose, a_transpose, A);

    // Compute the product: f'(q)*(f_x+g_y) + g'(q)*(f_x+g_y)
    for( int m1=1; m1 <= meqn; m1++ )
    {
        double tmp1 = 0.;
        double tmp2 = 0.;
        for( int m2=1; m2 <= meqn; m2++ )
        {
            tmp1 += -(A.get(1, m1, m2, 1)) * ( fx_val.get(m2) + gy_val.get(m2));
            tmp2 += -(A.get(1, m1, m2, 2)) * ( fx_val.get(m2) + gy_val.get(m2));
        }
        f_t.set( m1, tmp1 );
        g_t.set( m1, tmp2 );
    }

    // ---  Third-order terms --- //
    if( nterms > 2 )
    {

        // Hessian
        dTensor5 H( 1, meqn, meqn, meqn, 2 );
        D2FluxFunc(xpts, q_transpose, a_transpose, H);

        // ----------------------------------- //
        // Part I: Compute (f_x + g_y)_{,t}
        // ----------------------------------- //

        dTensor1 fx_plus_gy_t( meqn ); fx_plus_gy_t.setall(0.);

        // Start with a term that get's used frequently:
        dTensor1 fx_plus_gy( meqn );
        for( int m =1; m <= meqn; m++ )
        {
            double tmp = fx_val.get(m) + gy_val.get(m);
            fx_plus_gy.set(m, tmp );
        }

        // Second-derivatives:
        dTensor1 fxx_val( meqn ), gyy_val( meqn );
        Diff2( dx, Fvals,  fxx_val  );
        Diff2( dy, Gvals,  gyy_val  );

        // Cross - derivaties
        dTensor1 fxy_val  ( meqn );  fxy_val.setall(0.);
        dTensor1 gxy_val  ( meqn );  gxy_val.setall(0.);

        // Stencil for mixed derivatives
        //
        // TODO - this is a clunky way to compute these derivatives!
//      for( int m=1; m <= meqn; m++ )
//      {
//          dTensor1 tmpF(5);
//          dTensor1 tmpG(5);
//          for( int m1=-2; m1 <= 2; m1++ )
//          {
//              tmpF.set(m1+3, Diff1( dx, R.get(i-2,j+m1,m,1), R.get(i-1,j+m1,m,1), R.get(i,j+m1,m,1), R.get(i+1,j+m1,m,1), R.get(i+2,j+m1,m,1) ) );
//              tmpG.set(m1+3, Diff1( dx, R.get(i-2,j+m1,m,2), R.get(i-1,j+m1,m,2), R.get(i,j+m1,m,2), R.get(i+1,j+m1,m,2), R.get(i+2,j+m1,m,2) ) );
//          }
//          fxy_val.set(m, Diff1( dy, tmpF.get(1), tmpF.get(2), tmpF.get(3), tmpF.get(4), tmpF.get(5) ) );
//          gxy_val.set(m, Diff1( dy, tmpG.get(1), tmpG.get(2), tmpG.get(3), tmpG.get(4), tmpG.get(5) ) );
//      }

        // Clean, minimal stencil for computing u_xy using the smallest
        // fourth-order stencil available.
        for( int m=1; m <= meqn; m++ )
        {
            // Second-order terms
            double tmp = 0.25*(R.get(i+1,j+1,m,1) - R.get(i-1,j+1,m,1) - R.get(i+1,j-1,m,1) + R.get(i-1,j-1,m,1));
            // Higher-order terms:
            tmp -= (1./24.)*(
                R.get(i+2,j+1,m,1) + R.get(i-2,j-1,m,1) - R.get(i+2,j-1,m,1) - R.get(i-2,j+1,m,1) -
                R.get(i+1,j+2,m,1) - R.get(i-1,j-2,m,1) + R.get(i+1,j-2,m,1) + R.get(i-1,j+2,m,1) );
            tmp *= (1./(dx*dy));
            fxy_val.set(m, tmp );

            tmp = 0.25*(R.get(i+1,j+1,m,2) - R.get(i-1,j+1,m,2) - R.get(i+1,j-1,m,2) + R.get(i-1,j-1,m,2));
            tmp -= (1./24.)*(
                R.get(i+2,j+1,m,2) + R.get(i-2,j-1,m,2) - R.get(i+2,j-1,m,2) - R.get(i-2,j+1,m,2) -
                R.get(i+1,j+2,m,2) - R.get(i-1,j-2,m,2) + R.get(i+1,j-2,m,2) + R.get(i-1,j+2,m,2) );
            tmp *= (1./(dx*dy));
            gxy_val.set(m, tmp);

        }

        // Compute terms that get multiplied by 
        //     \pd2{ f }{ q } and \pd2{ g }{ q }.
        for( int m =1; m <= meqn; m++ )
        {
            double tmp = 0.;

            // Terms that get multiplied by the Hessian:
            for( int m1=1; m1 <= meqn; m1++ )
            for( int m2=1; m2 <= meqn; m2++ )
            {

                tmp += H.get(1,m,m1,m2,1)*qx_val.get(m1)*fx_plus_gy.get(m2);
                tmp += H.get(1,m,m1,m2,2)*qy_val.get(m1)*fx_plus_gy.get(m2);
            }

            // Terms that get multiplied by f'(q) and g'(q):
            for( int m1=1; m1 <= meqn; m1++ )
            {

                tmp += A.get(1,m,m1,1)*( fxx_val.get(m1)+gxy_val.get(m1) );
                tmp += A.get(1,m,m1,2)*( fxy_val.get(m1)+gyy_val.get(m1) );
            }

            fx_plus_gy_t.set( m, tmp );
        }


        // ----------------------------------- //
        // Part II: Compute 
        //      f'(q) * fx_plus_gy_t and 
        //      g'(q) * fx_plus_gy_t
        // ----------------------------------- //

        // Add in the third term that gets multiplied by A:
        for( int m1=1; m1 <= meqn; m1++ )
        {
            double tmp1 = 0.;
            double tmp2 = 0.;
            for( int m2=1; m2 <= meqn; m2++ )
            {
                tmp1 += A.get(1,m1,m2,1)*fx_plus_gy_t.get(m2);
                tmp2 += A.get(1,m1,m2,2)*fx_plus_gy_t.get(m2);
            }
            f_tt.set( m1, tmp1 );
            g_tt.set( m1, tmp2 );
        }

        // ----------------------------------------------- //
        // Part III: Add in contributions from
        //      f''(q) * (fx_plus_gy, fx_plus_gy ) and 
        //      g''(q) * (fx_plus_gy, fx_plus_gy ).
        // ----------------------------------------------- //
        for( int m =1; m <= meqn; m++ )
        {
            double tmp1 = 0.;
            double tmp2 = 0.;

            // Terms that get multiplied by the Hessian:
            for( int m1=1; m1 <= meqn; m1++ )
            for( int m2=1; m2 <= meqn; m2++ )
            {
                tmp1 += H.get(1,m,m1,m2,1)*fx_plus_gy.get(m1)*fx_plus_gy.get(m2);
                tmp2 += H.get(1,m,m1,m2,2)*fx_plus_gy.get(m1)*fx_plus_gy.get(m2);
            }

            f_tt.set( m, f_tt.get(m) + tmp1 );
            g_tt.set( m, g_tt.get(m) + tmp2 );
        }

    }
    else
    {
        // f_tt and g_tt not used.  Set to zero to be safe
        f_tt.setall(0.);
        g_tt.setall(0.);
    }

}

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1,
    const StateVars& Q1,
    double alpha2, double beta2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G)
{

    const dTensorBC3& q1   = Q1.const_ref_q  ();
    const dTensorBC3& aux1 = Q1.const_ref_aux();

    const dTensorBC3& q2   = Q2.const_ref_q  ();
    const dTensorBC3& aux2 = Q2.const_ref_aux();


    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double dy    = global_ini_params.get_dy();
    const double xlow  = global_ini_params.get_xlow();
    const double ylow  = global_ini_params.get_ylow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC4 R1( mx, my, meqn, 2, mbc );  // place-holder for the flux function
    dTensorBC4 R2( mx, my, meqn, 2, mbc );  // place-holder for the flux function
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q1, aux1, R1, &FluxFunc );
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q2, aux2, R2, &FluxFunc );

// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );

const int ndim = 2;

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    for( int j = 1-mbc_small; j <= my+mbc_small; j++ )
    {

        // Time derivatives of the flux function
        dTensor1 f1_t( meqn ), g1_t( meqn );
        dTensor1 f2_t( meqn ), g2_t( meqn );

        dTensor1 f1_tt( meqn ), g1_tt( meqn );
        dTensor1 f2_tt( meqn ), g2_tt( meqn );

        double xc = xlow + double(i)*dx - 0.5*dx;
        double yc = ylow + double(j)*dy - 0.5*dy;

        LocalIntegrate( 2, dx, dy, xc, yc, meqn, maux, mpts_sten, half_mpts_sten,
            i, j, q1, aux1, R1, f1_t, f1_tt, g1_t, g1_tt );

        LocalIntegrate( 2, dx, dy, xc, yc, meqn, maux, mpts_sten, half_mpts_sten,
            i, j, q2, aux2, R2, f2_t, f2_tt, g2_t, g2_tt );

        // Two-stage, two-derivative method:
        for( int m=1; m<=meqn; m++ )
        {
            F.set( i, j, m, 
                alpha1*R1.get(i,j,m,1) + dt*(beta1*f1_t.get(m)) + 
                alpha2*R2.get(i,j,m,1) + dt*(beta2*f2_t.get(m)) );

            G.set( i, j, m, 
                alpha1*R1.get(i,j,m,2) + dt*(beta1*g1_t.get(m)) +
                alpha2*R2.get(i,j,m,2) + dt*(beta2*g2_t.get(m)) );
        }

    }

}

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1, double charlie1, 
    const StateVars& Q1,
    double alpha2, double beta2, double charlie2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G)
{

    const dTensorBC3& q1   = Q1.const_ref_q  ();
    const dTensorBC3& aux1 = Q1.const_ref_aux();

    const dTensorBC3& q2   = Q2.const_ref_q  ();
    const dTensorBC3& aux2 = Q2.const_ref_aux();


    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double dy    = global_ini_params.get_dy();
    const double xlow  = global_ini_params.get_xlow();
    const double ylow  = global_ini_params.get_ylow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC4 R1( mx, my, meqn, 2, mbc );  // place-holder for the flux function
    dTensorBC4 R2( mx, my, meqn, 2, mbc );  // place-holder for the flux function
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q1, aux1, R1, &FluxFunc );
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q2, aux2, R2, &FluxFunc );

// TODO  - allow for different sized stencils for different orders (-DS)
const int mbc_small      = 3;
const int      mpts_sten = 5;
const int half_mpts_sten = (mbc+1)/2;    assert_eq( half_mpts_sten, 3 );

const int ndim = 2;

    // Compute finite difference approximations on all of the conserved
    // variables:
#pragma omp parallel for
    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    for( int j = 1-mbc_small; j <= my+mbc_small; j++ )
    {

        // Time derivatives of the flux function
        dTensor1 f1_t( meqn ), g1_t( meqn );
        dTensor1 f2_t( meqn ), g2_t( meqn );

        dTensor1 f1_tt( meqn ), g1_tt( meqn );
        dTensor1 f2_tt( meqn ), g2_tt( meqn );

        double xc = xlow + double(i)*dx - 0.5*dx;
        double yc = ylow + double(j)*dy - 0.5*dy;

        LocalIntegrate( 3, dx, dy, xc, yc, meqn, maux, mpts_sten, half_mpts_sten,
            i, j, q1, aux1, R1, f1_t, f1_tt, g1_t, g1_tt );

        LocalIntegrate( 3, dx, dy, xc, yc, meqn, maux, mpts_sten, half_mpts_sten,
            i, j, q2, aux2, R2, f2_t, f2_tt, g2_t, g2_tt );

        // Time-averaged flux function
        for( int m=1; m<=meqn; m++ )
        {

            F.set( i, j, m, 
                alpha1*R1.get(i,j,m,1) + dt*(beta1*f1_t.get(m) + charlie1*dt*f1_tt.get(m)) + 
                alpha2*R2.get(i,j,m,1) + dt*(beta2*f2_t.get(m) + charlie2*dt*f2_tt.get(m)) );

            G.set( i, j, m, 
                alpha1*R1.get(i,j,m,2) + dt*(beta1*g1_t.get(m) + charlie1*dt*g1_tt.get(m)) +
                alpha2*R2.get(i,j,m,2) + dt*(beta2*g2_t.get(m) + charlie2*dt*g2_tt.get(m)) );
        }

    }

}
#include <cmath>
#include "assert.h"            // for assert_eq.  Can be removed in future
#include "tensors.h"
#include "dog_math.h"
#include "StateVars.h"
#include "IniParams.h"
#include "ConstructL.h"

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x + g(q,x,t)_y = Psi(q,x,t)
//
void ConstructL( StateVars& Q, dTensorBC3& Lstar, dTensorBC3& smax)
{

    dTensorBC3&   q = Q.ref_q();
    dTensorBC3& aux = Q.ref_aux();

    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int    mbc = global_ini_params.get_mbc();

    // Size of the WENO stencil
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );

    // The flux, f_{i-1/2, j} and g_{i, j-1/2}.  Recall that the
    // flux lives at the nodal locations, i-1/2, so there is one more term in
    // that vector than on the original grid.  Additionally, for 2D code, there
    // are two terms in the flux: q_t + f_x + g_y = psi.
    dTensorBC3  Fhat(mx+1, my,   meqn, mbc );
    dTensorBC3  Ghat(mx,   my+1, meqn, mbc );

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();

    double alpha1 = 0.;
    double alpha2 = 0.;
    if( global_ini_params.get_global_alpha() )
    {
        // Global wave speed
        GlobalWaveSpd( Q, alpha1, alpha2);
    }

    // Normal vector.  This is a carry-over from the DG code.
    dTensor1 nvec(2);

    // --------------------------------------------------------------------- //
    // Compute Fhat{i-1/2, j} - 1st component of the flux function
    // --------------------------------------------------------------------- //
    nvec.set(1, 1.0 );  nvec.set(2, 0.0 );
#pragma omp parallel for
    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,j,m) + q.get(i-1,j,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(maux);
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i-1,j,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1         );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t(         ws+1, maux );

        dTensor2 xvals( ws+1, 2 );
        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int is = i-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( is )*dx - 0.5*dx;
            double yi = ylow + double( j  )*dy - 0.5*dy;
            xvals.set( s, 1, xi );
            xvals.set( s, 2, yi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(is, j, m ) );
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(is, j, ma ) );
            }
        }

        // The format of Flux and ProjectLeftEig/ProjectRightEig do not
        // contain the same order.  That is, Flux assumes q(1:npts, 1:meqn),
        // whereas the other functions assume q(1:meqn, 1:npts).  For
        // consistency, I will copy back to the latter, because the WENO
        // reconstruction *should* be faster if the list of points is second.
        // (-DS)
        ConvertTranspose( qvals,   qvals_t   );
        ConvertTranspose( auxvals, auxvals_t );

        // Sample the flux function over the stencil:
        dTensor3 fvals_t( ws+1, meqn, 2 );
        FluxFunc( xvals, qvals_t, auxvals_t, fvals_t );

        // Flux function "f" in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );
        for( int me=1; me <= meqn; me++ )
        for( int s=1; s <= ws+1; s++ )
        {
            f.set(me, s, fvals_t.get( s, me, 1 ) );  // 1st-component - f
            g.set(me, s, fvals_t.get( s, me, 2 ) );  // 2nd-component - g
        }


        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( 1, Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( 1, Auxavg, Qavg,     f, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------

        // -- Compute a local wave speed -- //

        dTensor1 xedge(2), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i-1, j, m) );
            Qr.set(m, q.get(i  , j, m) );
        }

        for( int m=1; m<= maux; m++)
        {
            Auxl.set(m, aux.get(i-1, j, m) );
            Auxr.set(m, aux.get(i  , j, m) );
        }

        // Compute an approximate "fastest" wave speed.
        // TODO - this is redundant in the case of a global value of alpha ...
        // (-DS 6/19/2014)
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( alpha1, Max( abs(s1), abs(s2) ) );
        smax.set( i, j, 1, Max( smax.get(i,j,1), alpha )  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here


        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
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
        ProjectRightEig(1, Auxavg, Qavg, ghat, fhat_loc);
        for( int m=1; m <= meqn; m++ )
        {
            Fhat.set(i, j, m, fhat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Compute Ghat{i, j-1/2} - 2nd-component of the flux function
    // --------------------------------------------------------------------- //
    nvec.set(1, 0.0 );  nvec.set(2, 1.0 );
#pragma omp parallel for
    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,j,m) + q.get(i,j-1,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(maux);
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i,j-1,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1         );
        dTensor2 qvals_t( ws+1, meqn ), auxvals_t(         ws+1, maux );
        dTensor2 xvals( ws+1, 2 );
        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int js = j-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( i  )*dx - 0.5*dx;
            double yi = ylow + double( js )*dy - 0.5*dy;
            xvals.set( s, 1, xi );
            xvals.set( s, 2, yi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(i, js, m ) );
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(i, js, ma ) );
            }
        }

        // The format of Flux and ProjectLeftEig/ProjectRightEig do not
        // contain the same order.  That is, Flux assumes q(1:npts, 1:meqn),
        // whereas the other functions assume q(1:meqn, 1:npts).  For
        // consistency, I will copy back to the latter, because the WENO
        // reconstruction *should* be faster if the list of points is second.
        // (-DS)
        ConvertTranspose( qvals,   qvals_t   );
        ConvertTranspose( auxvals, auxvals_t );

        // Sample f over the stencil:
        dTensor3 fvals_t( ws+1, meqn, 2 );
        FluxFunc( xvals, qvals_t, auxvals_t, fvals_t );

        // Flux function in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );
        for( int me=1; me <= meqn; me++ )
        for( int s=1; s <= ws+1; s++ )
        {
            f.set(me, s, fvals_t.get( s, me, 1 ) );
            g.set(me, s, fvals_t.get( s, me, 2 ) );
        }


        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( 2, Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( 2, Auxavg, Qavg,     g, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------

        // -- Compute a local wave speed -- //

        dTensor1 xedge(2), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i, j-1, m) );
            Qr.set(m, q.get(i, j,   m) );
        }

        for( int m=1; m<= maux; m++)
        {
            Auxl.set(m, aux.get(i, j-1, m) );
            Auxr.set(m, aux.get(i, j,   m) );
        }

        // Compute an approximate "fastest" wave speed.
        // TODO - this is redundant in the case of a global value of alpha ...
        // (-DS 6/19/2014)
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( alpha2, Max( abs(s1), abs(s2) ) );
        smax.set( i, j, 2, Max( smax.get(i,j,2), alpha )  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here

        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
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
        ProjectRightEig(2, Auxavg, Qavg, ghat, fhat_loc);
        for( int m=1; m <= meqn; m++ )
        {
            Ghat.set(i,j,m, fhat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_{i,j} = Lstar = -1/dx( fh_{i+1/2,j} - fh_{i-1/2,j} )
    //                           -1/dy( gh_{i,j+1/2} - gh_{i,j-1/2} ).
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {

        // Compute the source term.
        SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, Lstar, &SourceTermFunc);
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1, j,   m) - Fhat.get(i, j, m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,   j+1, m) - Ghat.get(i, j, m) ) / dy;
                Lstar.set(i,j, m, Lstar.get(i,j,m) + tmp );
            }
        }
    }
    else  // No source term
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1,j,  m) - Fhat.get(i,j,m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,  j+1,m) - Ghat.get(i,j,m) ) / dy;
                Lstar.set(i,j, m, tmp );
            }
        }
    }

    // ---------------------------------------------------------
    // Add extra contributions to Lstar
    // ---------------------------------------------------------
    // LstarExtra(node,aux,q,Lstar);

}

void ConvertTranspose( const dTensor2& qin, dTensor2& qout )
{
    const int m1 = qin.getsize(1);
    const int m2 = qin.getsize(2);
    assert_eq( m1, qout.getsize(2) );
    assert_eq( m2, qout.getsize(1) );

    for( int i=1; i<= m1; i++ )
    for( int j=1; j<= m2; j++ )
    {
        qout.set(j,i, qin.get(i,j) );
    }

}
#ifndef _CONSTRUCT_L_H_
#define _CONSTRUCT_L_H_

#include "StateVars.h"
#include "tensors.h"

// --- User supplied functions --- //
void FluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor3& flux);
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, 
            const dTensor2& auxvals, dTensor2& source);

void ProjectLeftEig( int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
    const dTensor2& Qvals, dTensor2& Wvals);
void ProjectRightEig(int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
                     const dTensor2& Wvals, dTensor2& Qvals);

// Routines for computing maximum wave speeds

// Local wave speed
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge, 
    const dTensor1& Ql,   const dTensor1& Qr, 
    const dTensor1& Auxl, const dTensor1& Auxr,
    double& s1,double& s2);

// Global wave speed
void GlobalWaveSpd( const StateVars& Q, double& alpha1, double& alpha2);

// Routine for sampling a given function
void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  dTensorBC3& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

// Routine to deal with the silly mess where the Fluxes and the
// Projections are all defined separately.
void ConvertTranspose( const dTensor2& qin, dTensor2& qout );

#endif
#include <iostream>
#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"
#include "assert.h"
#include "StateVars.h"

using namespace std;

// User supplied functions defining the Flux function, Jacobian, and
// Hessian of the flux function.
void FluxFunc(const dTensor2&,const dTensor2&,const dTensor2&,dTensor3&);

// Used for construcing the flux function
void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC4& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&));

// Same as above, but this one only returns flux values, and not just the RHS.
void ConstructLFL( const double dt, StateVars& Qnew, 
    dTensorBC3& fLF, dTensorBC3& gLF, dTensorBC3& Lstar, dTensorBC3& smax )
{

    dTensorBC3&    q = Qnew.ref_q  ();
    dTensorBC3&  aux = Qnew.ref_aux();

    void SetBndValues(StateVars& Qnew); 
    SetBndValues( Qnew );

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC4 R( mx, my, meqn, 2, mbc );  // place-holder for the flux function

    const int mbc_small      = 3;

    // Compute finite difference approximations on all of the conserved
    // variables:
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, R, &FluxFunc );

    // Find a global wave speed
    void GlobalWaveSpd( const StateVars& Q, double& alpha1, double& alpha2);
    double alpha_x, alpha_y;
    GlobalWaveSpd( Qnew, alpha_x, alpha_y );
    smax.set(1,1,1, alpha_x);
    smax.set(1,1,2, alpha_y);
//  printf("alpha_x, alpha_y = %f, %f \n", alpha_x, alpha_y );

/* This section of code was written specifically for Euler equations */
/*
    alpha_x = 1.0e-15;
    alpha_y = 1.0e-15;
    const double gamma = global_ini_params.get_gamma();

    for (int i = 1; i <= mx;  i++)
    for (int j = 1; j <= my;  j++)
    {
        // -- Compute a local wave speed //
        dTensor1 xedge(2);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        const double rho = q.get(i,j,1);
        const double u1 =  q.get(i,j,2)/rho;
        const double u2 =  q.get(i,j,3)/rho;
        const double u3 =  q.get(i,j,4)/rho;
        const double energy = q.get(i,j,5);

        const double press = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

        if (press < 0.0)
            cout << " x negative press = " << press << " (x,y)= "<< 
            xedge.get(1)<<" "<<xedge.get(2) << endl;

        const double c = sqrt(fabs(gamma*press/rho));
        const double alpha_xl = abs(u1) + c;
        const double alpha_yl = abs(u2) + c;

        smax.set( i, j, 1, Max( smax.get(i,j,1), alpha_xl )  );
        smax.set( i, j, 2, Max( smax.get(i,j,2), alpha_yl )  );

        alpha_x = Max(alpha_x, alpha_xl);
        alpha_y = Max(alpha_y, alpha_yl);
    }
*/

    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =  0.5*(R.get(i-1,j,m,1)+R.get(i,j,m,1)) 
                   + 0.5*alpha_x*(q.get(i-1,j,m)-q.get(i,j,m));
        fLF.set(i,j,m, hf );
    }

    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =  0.5*(R.get(i,j-1,m,2)+R.get(i,j,m,2)) 
                   + 0.5*alpha_y*(q.get(i,j-1,m)-q.get(i,j,m));
        gLF.set(i,j,m, hf );
    }

    // Construct Lstar, defined by:
    //
    //    d/dt q_{i,j} = Lstar = -1/dx( fh_{i+1/2,j} - fh_{i-1/2,j} )
    //                           -1/dy( gh_{i,j+1/2} - gh_{i,j-1/2} ).
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {
        printf("Error: source-term not implemented for Lax-Wendroff method\n");
        exit(1);
    }
    else  // No source term
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(fLF.get(i+1,j,  m) - fLF.get(i,j,m) ) / dx;
                tmp =  tmp   -(gLF.get(i,  j+1,m) - gLF.get(i,j,m) ) / dy;
                Lstar.set(i,j,m, tmp );
            }
        }
    }

}

// Construct a Lax-Friedrich's flux
//
// This is essentially a clone of the above function.  The difference here is
// that the boundary conditions are applied twice, once before looping over F,
// and once before looping over g.  This is necessary for the problems with
// "geometry".
//
void ConstructLFL_BCx2( const double dt, StateVars& Qnew,
    dTensorBC3& Fhat, dTensorBC3& Ghat,
    dTensorBC3& Lstar, dTensorBC3& smax)

{

    dTensorBC3&    q = Qnew.ref_q  ();
    dTensorBC3&  aux = Qnew.ref_aux();

    void SetBndValuesX(StateVars& Q); // Only set conditions along x-direction
    void SetBndValuesY(StateVars& Q); // Only set conditions along y-direction

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();

    const int mbc_small      = 3;

    // Sample the flux function on the entire domain:
    dTensorBC4 R( mx, my, meqn, 2, mbc );  // place-holder for the flux function

    // Compute finite difference approximations on all of the conserved
    // variables:
    SetBndValuesX(Qnew);
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, R, &FluxFunc );

    // Find a global wave speed
    void GlobalWaveSpd( const StateVars& Q, double& alpha1, double& alpha2);
    double alpha_x, alpha_y;
    GlobalWaveSpd( Qnew, alpha_x, alpha_y );

    // Find fastest wave speed in X-direction
    //
    // This chunk of code was written specifically for Euler equations
    //
//  double alpha_x     = 1.0e-15;
//  const double gamma = global_ini_params.get_gamma();
//  for (int i = 1; i <= mx;  i++)
//  for (int j = 1; j <= my;  j++)
//  {
//      // -- Compute a local wave speed -x- //
//      dTensor1 xedge(2);
//      xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
//      xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

//      const double rho = q.get(i,j,1);
//      const double u1 =  q.get(i,j,2)/rho;
//      const double u2 =  q.get(i,j,3)/rho;
//      const double u3 =  q.get(i,j,4)/rho;
//      const double energy = q.get(i,j,5);
//      const double press = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

//      if (press < 0.0)
//      {
//          cout << " x negative press = " << press << " (x,y)= "<< 
//          xedge.get(1)<<" "<<xedge.get(2) << endl;
//      }

//      const double c = sqrt(fabs(gamma*press/rho));
//      const double alpha = abs(u1) + c;

//      smax.set( i, j, 1, Max( smax.get(i,j,1), alpha )  );
//      alpha_x = Max(alpha_x, alpha);
//  }

    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =  0.5*(R.get(i-1,j,m,1)+R.get(i,j,m,1)) 
                   + 0.5*alpha_x*(q.get(i-1,j,m)-q.get(i,j,m));
        Fhat.set(i,j,m, hf );
    }

//  SetBndValuesY(Qnew);
//  SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, R, &FluxFunc );

    // Find fastest wave speed in Y-direction
    //
    // This chunk of code was written specifically for Euler equations
    //
//  double alpha_y = 1.0e-15;
//  for (int i = 1; i <= mx;  i++)
//  for (int j = 1; j <= my;  j++)
//  {
//      // -- Compute a local wave speed -y- //
//      dTensor1 xedge(2); 
//      xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
//      xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

//      const double rho = q.get(i,j,1);
//      const double u1 =  q.get(i,j,2)/rho;
//      const double u2 =  q.get(i,j,3)/rho;
//      const double u3 =  q.get(i,j,4)/rho;
//      const double energy = q.get(i,j,5);

//      const double press = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

//      if (press < 0.0)
//      {
//          cout << " y negative press = " << press;
//          cout << " (x,y)= "<< xedge.get(1)<<" "<<xedge.get(2) << endl;
//      }

//      // Sound speeds
//      const double c = sqrt(fabs(gamma*press/rho));
//      const double alpha = abs(u2) + c;

//      smax.set( i, j, 2, Max( smax.get(i,j,2), alpha )  );
//      alpha_y = Max(alpha_y, alpha);
//  }

#pragma omp parallel for
    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =   0.5*(R.get(i,j-1,m,2)+R.get(i,j,m,2)) 
                    + 0.5*alpha_y*(q.get(i,j-1,m)-q.get(i,j,m));
        Ghat.set(i,j,m, hf );
    }

    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_{i,j} = Lstar = -1/dx( fh_{i+1/2,j} - fh_{i-1/2,j} )
    //                           -1/dy( gh_{i,j+1/2} - gh_{i,j-1/2} ).
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {
        printf("Error: source-term not implemented for Lax-Wendroff method\n");
        exit(1);
    }
    else  // No source term
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1,j,  m) - Fhat.get(i,j,m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,  j+1,m) - Ghat.get(i,j,m) ) / dy;
                Lstar.set(i,j,m, tmp );
            }
        }
    }

}

#include <cmath>
#include "assert.h"            // for assert_eq.  Can be removed in future
#include "tensors.h"
#include "dog_math.h"
#include "IniParams.h"
#include "StateVars.h"

// --- User supplied functions --- //
void ProjectLeftEig( int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
    const dTensor2& Qvals, dTensor2& Wvals);
void ProjectRightEig(int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
                     const dTensor2& Wvals, dTensor2& Qvals);
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge, 
    const dTensor1& Ql,   const dTensor1& Qr, 
    const dTensor1& Auxl, const dTensor1& Auxr,
    double& s1,double& s2);
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, 
            const dTensor2& auxvals, dTensor2& source);

void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  dTensorBC3& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));


// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x + g(q,x,t)_y = Psi(q,x,t)
//
// EXPERIMENTAL CODE - This routine performs the Lax-Friedrich's flux
// splitting on a modified flux function, F and G.
void ConstructLxWL( const StateVars& Q,
        dTensorBC3& F,         // <--- new term: integrated flux, f
        dTensorBC3& G,         // <--- new term: integrated flux, g
        dTensorBC3& Lstar,
        dTensorBC3& smax)
{

    const dTensorBC3&    q = Q.const_ref_q  ();
    const dTensorBC3&  aux = Q.const_ref_aux();
    

    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Routine to deal with the silly mess where the Fluxes and the
    // Projections are all defined separately.
    void ConvertTranspose( const dTensor2& qin, dTensor2& qout );

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int    mbc = global_ini_params.get_mbc();

    // Size of the WENO stencil
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );

    // The flux, f_{i-1/2, j} and g_{i, j-1/2}.  Recall that the
    // flux lives at the nodal locations, i-1/2, so there is one more term in
    // that vector than on the original grid.  Additionally, for 2D code, there
    // are two terms in the flux: q_t + f_x + g_y = psi.
    dTensorBC3  Fhat(mx+1, my,   meqn, mbc );
    dTensorBC3  Ghat(mx,   my+1, meqn, mbc );

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();

    // Terms used in the case of using a global alpha
    double alpha1 = 0.;
    double alpha2 = 0.;
    if( global_ini_params.get_global_alpha() )
    {
        // Global wave speed
        void GlobalWaveSpd( const StateVars& Q, 
            double& alpha1, double& alpha2);
        GlobalWaveSpd( Q, alpha1, alpha2);
    }

    // Normal vector.  This is a carry-over from the DG code.
    dTensor1 nvec(2);

    // --------------------------------------------------------------------- //
    // Compute F{i-1/2, j} - 1st component of the flux function
    // --------------------------------------------------------------------- //
    nvec.set(1, 1.0 );  nvec.set(2, 0.0 );
#pragma omp parallel for
    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,j,m) + q.get(i-1,j,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg( maux );
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i-1,j,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1         );

        // Flux function "f" and "g" in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );
        dTensor2 xvals( ws+1, 2 );

        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int is = i-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( is )*dx - 0.5*dx;
            double yi = ylow + double( j  )*dy - 0.5*dy;
            xvals.set( s, 1, xi );
            xvals.set( s, 2, yi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(is, j, m ) );
                f.set    ( m, s, F.get(is, j, m ) );  // <-- NEW part (sample integrated flux)
                g.set    ( m, s, G.get(is, j, m ) );  // <-- NEW part (sample integrated flux)
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(is, j, ma ) );
            }
        }

        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( 1, Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( 1, Auxavg, Qavg,     f, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------

        // -- Compute a local wave speed -- //

        dTensor1 xedge(2), Ql(meqn), Qr(meqn);
        dTensor1 Auxl( maux ), Auxr( maux );
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i-1, j, m) );
            Qr.set(m, q.get(i  , j, m) );
        }

        for( int m=1; m<=maux; m++)
        {
            Auxl.set(m, aux.get(i-1, j, m) );
            Auxr.set(m, aux.get(i  , j, m) );
        }

        // Compute an approximate "fastest" wave speed.
        // TODO - this is redundant in the case of a global value of alpha ...
        // (-DS 6/19/2014)
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( alpha1, Max( abs(s1), abs(s2) ) );
        smax.set( i, j, 1, Max( smax.get(i,j,1), alpha )  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here


        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
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
        ProjectRightEig(1, Auxavg, Qavg, ghat, fhat_loc);
        for( int m=1; m <= meqn; m++ )
        {
            Fhat.set(i, j, m, fhat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Compute G{i, j-1/2} - 2nd-component of the flux function
    // --------------------------------------------------------------------- //
    nvec.set(1, 0.0 );  nvec.set(2, 1.0 );
#pragma omp parallel for
    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,j,m) + q.get(i,j-1,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(maux );
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i,j-1,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1         );
        dTensor2 xvals( ws+1, 2 );

        // Flux function in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );

        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int js = j-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( i  )*dx - 0.5*dx;
            double yi = ylow + double( js )*dy - 0.5*dy;
            xvals.set( s, 1, xi );
            xvals.set( s, 2, yi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(i, js, m ) );
                f.set(m,s, F.get(i,js,m ) );  // <-- NEW part (sample integrated flux)
                g.set(m,s, G.get(i,js,m ) );  // <-- NEW part (sample integrated flux)
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(i, js, ma ) );
            }
        }

        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( 2, Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( 2, Auxavg, Qavg,     g, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------

        // -- Compute a local wave speed -- //

        dTensor1 xedge(2), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );
        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i, j-1, m) );
            Qr.set(m, q.get(i, j,   m) );
        }

        for( int m=1; m<= maux; m++)
        {
            Auxl.set(m, aux.get(i, j-1, m) );
            Auxr.set(m, aux.get(i, j,   m) );
        }

        // Compute an approximate "fastest" wave speed.
        // TODO - this is redundant in the case of a global value of alpha ...
        // (-DS 6/19/2014)
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( alpha2, Max( abs(s1), abs(s2) ) );
        smax.set( i, j, 2, Max( smax.get(i,j,2), alpha )  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here

        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
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

        dTensor2 ghat_loc( ghat );
        ProjectRightEig(2, Auxavg, Qavg, ghat, ghat_loc );
        for( int m=1; m <= meqn; m++ )
        {
            Ghat.set(i,j,m, ghat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_{i,j} = Lstar = -1/dx( fh_{i+1/2,j} - fh_{i-1/2,j} )
    //                           -1/dy( gh_{i,j+1/2} - gh_{i,j-1/2} ).
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {
        printf("Error: source-term not implemented for Lax-Wendroff method\n");
        exit(1);
//      // Compute the source term.
//      SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, Lstar, &SourceTermFunc);
//#pragma omp parallel for
//      for (int i=1; i<=mx; i++)
//      for (int j=1; j<=my; j++)
//      {
//          for (int m=1; m<=meqn; m++)
//          {
//              double tmp = -(Fhat.get(i+1, j,   m) - Fhat.get(i, j, m) ) / dx;
//              tmp =  tmp   -(Ghat.get(i,   j+1, m) - Ghat.get(i, j, m) ) / dy;
//              Lstar.set(i,j, m, Lstar.get(i,j,m) + tmp );
//          }
//      }
    }
    else  // No source term
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1,j,  m) - Fhat.get(i,j,m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,  j+1,m) - Ghat.get(i,j,m) ) / dy;
                Lstar.set(i,j,m, tmp );
            }
        }
    }

    // ---------------------------------------------------------
    // Add extra contributions to Lstar
    // ---------------------------------------------------------
    // LstarExtra(aux,q,Lstar);

}

// Clone of the above function - however, this one saves flux values
void ConstructLxWL( const StateVars& Q,
        dTensorBC3& F, dTensorBC3& G,
        dTensorBC3& Fhat, dTensorBC3& Ghat,
        dTensorBC3& Lstar, dTensorBC3& smax)
{

    const dTensorBC3&    q = Q.const_ref_q  ();
    const dTensorBC3&  aux = Q.const_ref_aux();
    

    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Routine to deal with the silly mess where the Fluxes and the
    // Projections are all defined separately.
    void ConvertTranspose( const dTensor2& qin, dTensor2& qout );

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int    mbc = global_ini_params.get_mbc();

    // Size of the WENO stencil
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();

    // Terms used in the case of using a global alpha
    double alpha1 = 0.;
    double alpha2 = 0.;
    if( global_ini_params.get_global_alpha() )
    {
        // Global wave speed
        void GlobalWaveSpd( const StateVars& Q, 
            double& alpha1, double& alpha2);
        GlobalWaveSpd( Q, alpha1, alpha2);
    }


    // Normal vector.  This is a carry-over from the DG code.
    dTensor1 nvec(2);

    // --------------------------------------------------------------------- //
    // Compute F{i-1/2, j} - 1st component of the flux function
    // --------------------------------------------------------------------- //
    nvec.set(1, 1.0 );  nvec.set(2, 0.0 );
#pragma omp parallel for
    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,j,m) + q.get(i-1,j,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg( maux );
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i-1,j,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1         );

        // Flux function "f" and "g" in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );
        dTensor2 xvals( ws+1, 2 );

        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int is = i-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( is )*dx - 0.5*dx;
            double yi = ylow + double( j  )*dy - 0.5*dy;
            xvals.set( s, 1, xi );
            xvals.set( s, 2, yi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(is, j, m ) );
                f.set    ( m, s, F.get(is, j, m ) );  // <-- NEW part (sample integrated flux)
                g.set    ( m, s, G.get(is, j, m ) );  // <-- NEW part (sample integrated flux)
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(is, j, ma ) );
            }
        }

        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( 1, Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( 1, Auxavg, Qavg,     f, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------

        // -- Compute a local wave speed -- //

        dTensor1 xedge(2), Ql(meqn), Qr(meqn);
        dTensor1 Auxl( maux ), Auxr( maux );
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i-1, j, m) );
            Qr.set(m, q.get(i  , j, m) );
        }

        for( int m=1; m<=maux; m++)
        {
            Auxl.set(m, aux.get(i-1, j, m) );
            Auxr.set(m, aux.get(i  , j, m) );
        }

        // Compute an approximate "fastest" wave speed.
        // TODO - this is redundant in the case of a global value of alpha ...
        // (-DS 6/19/2014)
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( alpha1, Max( abs(s1), abs(s2) ) );
        smax.set( i, j, 1, Max( smax.get(i,j,1), alpha )  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here


        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
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
        ProjectRightEig(1, Auxavg, Qavg, ghat, fhat_loc);
        for( int m=1; m <= meqn; m++ )
        {
            Fhat.set(i, j, m, fhat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Compute G{i, j-1/2} - 2nd-component of the flux function
    // --------------------------------------------------------------------- //
    nvec.set(1, 0.0 );  nvec.set(2, 1.0 );
#pragma omp parallel for
    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
    {

        // --------------------------------------------------------------------
        // Part I: Compute Roe Averages
        //         TODO - the User may want to replace this ...
        // --------------------------------------------------------------------
        dTensor1 Qavg(meqn);
        for( int m=1; m <= meqn; m++ )
        {
            double tmp = 0.5*( q.get(i,j,m) + q.get(i,j-1,m) );
            Qavg.set(m, tmp );
        }
        dTensor1 Auxavg(maux );
        for( int ma=1; ma <= maux; ma++ )
        {
            double tmp = 0.5*( aux.get(i,j,ma) + aux.get(i,j-1,ma) );
            Auxavg.set(ma, tmp );
        }

        // --------------------------------------------------------------------
        // Part II: Compute w_{i+r} = R^{-1} q  and g_{i+r} = R^{-1} f
        // --------------------------------------------------------------------

        // Sample q over the stencil:
        dTensor2  qvals( meqn, ws+1  ), auxvals  ( maux, ws+1         );
        dTensor2 xvals( ws+1, 2 );

        // Flux function in q_t + f_x + g_y = 0:
        dTensor2 f( meqn, ws+1 ), g( meqn, ws+1 );

        for( int s=1; s <= ws+1; s++ )
        {
            // Index into the large array
            int js = j-1+s-r;

            // TODO - check that this is the correct value ...
            double xi = xlow + double( i  )*dx - 0.5*dx;
            double yi = ylow + double( js )*dy - 0.5*dy;
            xvals.set( s, 1, xi );
            xvals.set( s, 2, yi );

            for( int m=1; m <= meqn; m++ )
            {
                qvals.set( m, s, q.get(i, js, m ) );
                f.set(m,s, F.get(i,js,m ) );  // <-- NEW part (sample integrated flux)
                g.set(m,s, G.get(i,js,m ) );  // <-- NEW part (sample integrated flux)
            }
            for( int ma=1; ma <= maux; ma++ )
            {
                auxvals.set( ma, s, aux.get(i, js, ma ) );
            }
        }

        // Project entire stencil onto the characteristic variables:
        dTensor2 wvals( meqn, ws+1  ), gvals( meqn, ws+1 );
        ProjectLeftEig( 2, Auxavg, Qavg, qvals, wvals );
        ProjectLeftEig( 2, Auxavg, Qavg,     g, gvals );

        // --------------------------------------------------------------------
        // Part III: Apply Lax-Friedrich's flux splitting to g
        // --------------------------------------------------------------------

        // -- Compute a local wave speed -- //

        dTensor1 xedge(2), Ql(meqn), Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );
        for( int m=1; m<= meqn; m++)
        {
            Ql.set(m, q.get(i, j-1, m) );
            Qr.set(m, q.get(i, j,   m) );
        }

        for( int m=1; m<= maux; m++)
        {
            Auxl.set(m, aux.get(i, j-1, m) );
            Auxr.set(m, aux.get(i, j,   m) );
        }

        // Compute an approximate "fastest" wave speed.
        // TODO - this is redundant in the case of a global value of alpha ...
        // (-DS 6/19/2014)
        double s1,s2;
        SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

        const double alpha = Max( alpha2, Max( abs(s1), abs(s2) ) );
        smax.set( i, j, 2, Max( smax.get(i,j,1), alpha )  );
        const double l_alpha = global_ini_params.get_alpha_scaling()*alpha;  // extra safety factor added here

        // -- Flux splitting -- //

        dTensor2 gp( meqn, ws ), gm( meqn, ws );
        for( int m=1; m <= meqn; m++ )
        for( int s=1; s <= ws; s++ )
        {
            gp.set( m, s, 0.5*(gvals.get(m,s     ) + l_alpha*wvals.get(m,s     ) ) );
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

        dTensor2 ghat_loc( ghat );
        ProjectRightEig(2, Auxavg, Qavg, ghat, ghat_loc );
        for( int m=1; m <= meqn; m++ )
        {
            Ghat.set(i,j,m, ghat_loc.get(m,1) );
        }

    }
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_{i,j} = Lstar = -1/dx( fh_{i+1/2,j} - fh_{i-1/2,j} )
    //                           -1/dy( gh_{i,j+1/2} - gh_{i,j-1/2} ).
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {
        printf("Error: source-term not implemented for Lax-Wendroff method\n");
        exit(1);
//      // Compute the source term.
//      SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, Lstar, &SourceTermFunc);
//#pragma omp parallel for
//      for (int i=1; i<=mx; i++)
//      for (int j=1; j<=my; j++)
//      {
//          for (int m=1; m<=meqn; m++)
//          {
//              double tmp = -(Fhat.get(i+1, j,   m) - Fhat.get(i, j, m) ) / dx;
//              tmp =  tmp   -(Ghat.get(i,   j+1, m) - Ghat.get(i, j, m) ) / dy;
//              Lstar.set(i,j, m, Lstar.get(i,j,m) + tmp );
//          }
//      }
    }
    else  // No source term
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1,j,  m) - Fhat.get(i,j,m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,  j+1,m) - Ghat.get(i,j,m) ) / dy;
                Lstar.set(i,j,m, tmp );
            }
        }
    }

    // ---------------------------------------------------------
    // Add extra contributions to Lstar
    // ---------------------------------------------------------
    // LstarExtra(aux,q,Lstar);

}
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "IniParams.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "hooks.h"              // hooks (files that a user may wish to relink)
#include "constructs.h"
#include "app_defined.h"        // application (required) files
#include "misc2d.h"
#include "StateVars.h"
#include "FinSolveLxW.h"

using namespace std;

void FinSolveLxW( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC3& qnew = Qnew.ref_q  ();
    dTensorBC3&  aux = Qnew.ref_aux();

    // Time stepping information
    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number
    double t                  = Qnew.get_t();
    double dt                 = dtv[1];   // Start with time step from last frame
    double cfl                = 0.0;      // current CFL number
    double dtmin              = dt;       // Counters for max and min time step taken
    double dtmax              = dt;

    // Grid information
    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();
    const int meqn = global_ini_params.get_meqn();
    const int maux = global_ini_params.get_maux();
    const int mbc  = global_ini_params.get_mbc();
    const int numel = qnew.numel();

    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

    // Maximum wave speed at each flux interface value.  Note that the size of
    // this tensor is one longer in each direction.
    dTensorBC3 smax( mx+1, my+1, 2, mbc );           

    // Needed for rejecting a time step
    StateVars Qold( t, mx, my, meqn, maux, mbc );
    dTensorBC3& qold   = Qold.ref_q();
    dTensorBC3& auxold = Qold.ref_aux();
    Qold.copyfrom( Qnew );

    // Intermediate stages
    StateVars Qstar( t, mx, my, meqn, maux, mbc );
    dTensorBC3&   qstar = Qstar.ref_q();
    dTensorBC3& auxstar = Qstar.ref_aux();
    Qstar.copyfrom( Qnew );

    // Allocate storage for this solver
    dTensorBC3       F(mx, my, meqn, mbc );  // time-integrated flux
    dTensorBC3       G(mx, my, meqn, mbc );  // time-integrated flux
    dTensorBC3   Lstar(mx, my, meqn, mbc);   // Right hand side of ODE

    // Storage for the MPP limiter
    dTensorBC3* fhat;
    dTensorBC3* fLF;
    dTensorBC3* ghat;
    dTensorBC3* gLF;
    if( global_ini_params.get_mpp_limiter() )
    {
        fhat = new dTensorBC3( mx+1, my, meqn, mbc );
        fLF  = new dTensorBC3( mx+1, my, meqn, mbc );

        ghat = new dTensorBC3( mx, my+1, meqn, mbc );
        gLF  = new dTensorBC3( mx, my+1, meqn, mbc );
    }


    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;                           // Number of time steps taken
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step > nv )
        {
            printf(" Error in FinSolveRK.cpp: "         );
            printf("Exceeded allowed # of time steps \n");
            printf("    n_step = %d\n", n_step          );
            printf("        nv = %d\n", nv              );
            printf("Terminating program.\n"             );
            exit(1);
        }        

        // copy qnew into qold
        Qold.copyfrom( Qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while( m_accept==0 )
        {

            // set current time
            Qnew.set_t( t );
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, Qold, Qnew);

            // ---------------------------------------------------------
            // Take a full time step of size dt
            BeforeStep(dt, Qnew);
            SetBndValues(Qnew);
            ConstructIntegratedR( dt, Qnew, smax, F, G);


            if( global_ini_params.get_mpp_limiter() )
            {

                // Construct the high-order flux
                ConstructLxWL( Qnew, F, G, *fhat, *ghat, Lstar, smax );
       
                // Construct the low-order flux
                ConstructLFL( dt, Qnew, *fLF, *gLF, Lstar, smax );

                // Limit the high-order flux
                ApplyMPPLimiter2D( dt, qnew, *fLF, *gLF, *fhat, *ghat );

                // Update the solution:
#pragma omp parallel for
                for( int i=1; i <= mx; i++   )
                for( int j=1; j <= my; j++   )
                for( int m=1; m <= meqn; m++ )
                {
                    double tmp = (fhat->get(i+1,j,m)-fhat->get(i,j,m) );
                    qnew.set(i, j, m, qnew.get(i,j,m) - (dt/dx)*tmp );
                    tmp = (ghat->get(i,j+1,m)-ghat->get(i,j,m) );
                    qnew.set(i, j, m, qnew.get(i,j,m) - (dt/dy)*tmp );

// Test the LF solver by uncommenting this chunk
//                  double tmp = Lstar.get(i,j,m);
//                  qnew.set(i, j, m, qnew.get(i,j,m) + dt*tmp );

                }

            }
            else
            { 

                // Construct RHS
                ConstructLxWL( Qnew, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }

            }

            Qnew.set_t( Qnew.get_t() + dt );

            // Perform any extra work required:
            AfterStep(dt, Qnew);
            // ---------------------------------------------------------

            // do any extra work
            AfterFullTimeStep(dt, Qold, Qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveLxW2D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)       // accept
            { m_accept = 1; }
            else                    //reject
            {   
                t = told;
                if( global_ini_params.get_verbosity() )
                {
                    cout<<"FinSolveLxW2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                Qnew.copyfrom( Qold  );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln(Qnew);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    // Clean up allocated memory
    if( global_ini_params.get_mpp_limiter() )
    {
        delete fhat;
        delete fLF;
        delete ghat;
        delete gLF;
    }

}
#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

#include "StateVars.h"

// ------------------------------------------------------------
// Function definitions
void ConSoln( const StateVars& Q );

double GetCFL(double dt, double dtmax, const dTensorBC3& aux, const dTensorBC3& smax);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, StateVars& Q );
void SetBndValues(StateVars& Q );
void AfterStep(double dt, StateVars& Q);

// Called once before each full time step (c.f. BeforeStep, which is called
// before each stage in an RK method)
void BeforeFullTimeStep(double dt, StateVars& Qold, StateVars& Q );
void AfterFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew );

// ------------------------------------------------------------
// Lax-Wendroff information:
//
// TODO - copy over the ConstructIntegratedF over to this 
//        routine
//
// ------------------------------------------------------------

void ConstructLxWL( const StateVars& Q,
        dTensorBC3& F, dTensorBC3& G,
        dTensorBC3& Lstar, dTensorBC3& smax);

void ConstructLxWL( const StateVars& Q,
        dTensorBC3& F, dTensorBC3& G,
        dTensorBC3& Fhat, dTensorBC3& Ghat,
        dTensorBC3& Lstar, dTensorBC3& smax);

void ConstructIntegratedR( double dt, const StateVars& Q, 
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

// Construct a Lax-Friedrich's flux
void ConstructLFL( const double dt, StateVars& Q,
    dTensorBC3& Fhat, dTensorBC3& Ghat,
    dTensorBC3& Lstar, dTensorBC3& smax);

void ApplyMPPLimiter2D( 
        const double dt, const dTensorBC3& q, 
        const dTensorBC3& fLF, const dTensorBC3& gLF,
        dTensorBC3& fHat, dTensorBC3& gHat );

#endif
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "FinSolveLxW.h"     // functions directly called from this function
#include "constructs.h" 
#include "StateVars.h"
#include "IniParams.h"

// ------------------------------------------------------------
// Multiderivative integration
//
// These functions are for the two-stage methods.  One contains
// two-derivatives, and the second contains three derivatives.
void ConstructIntegratedR( double dt, 
    double alpha1, double beta1,
    const StateVars& Q1,
    double alpha2, double beta2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1, double charlie1,
    const StateVars& Q1,
    double alpha2, double beta2, double charlie2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);
// ------------------------------------------------------------


using namespace std;

void FinSolveMD( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC3& qnew = Qnew.ref_q  ();
    dTensorBC3&  aux = Qnew.ref_aux();

    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();

    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number
    double t                  = Qnew.get_t();
    double dt                 = dtv[1];   // Start with time step from last frame
    double cfl                = 0.0;      // current CFL number
    double dtmin              = dt;       // Counters for max and min time step taken
    double dtmax              = dt;

    // Grid information
    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();
    const int meqn = global_ini_params.get_meqn();
    const int maux = global_ini_params.get_maux();
    const int mbc  = global_ini_params.get_mbc();
    const int numel = qnew.numel();

    // Maximum wave speed at each flux interface value.  Note that the size of
    // this tensor is one longer in each direction.
    dTensorBC3 smax( mx+1, my+1, 2, mbc );

    // Needed for rejecting a time step
    StateVars Qold( t, mx, my, meqn, maux, mbc );
    dTensorBC3& qold   = Qold.ref_q();
    dTensorBC3& auxold = Qold.ref_aux();
    Qold.copyfrom( Qnew );

    // Intermediate stages
    StateVars Qstar( t, mx, my, meqn, maux, mbc );
    dTensorBC3&   qstar = Qstar.ref_q();
    dTensorBC3& auxstar = Qstar.ref_aux();
    Qstar.copyfrom( Qnew );


    // Allocate storage for this solver
    dTensorBC3   Lstar(mx, my, meqn, mbc);   // right hand side (for Euler steps)
    dTensorBC3       F(mx, my, meqn, mbc );  // time-integrated flux
    dTensorBC3       G(mx, my, meqn, mbc );  // time-integrated flux

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step > nv )
        {
            printf(" Error in FinSolveRK.cpp: "         );
            printf("Exceeded allowed # of time steps \n");
            printf("    n_step = %d\n", n_step          );
            printf("        nv = %d\n", nv              );
            printf("Terminating program.\n"             );
            exit(1);
        }        

        // copy qnew into qold
        Qold.copyfrom( Qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while( m_accept==0 )
        {

            // set current time
            Qnew.set_t( t );
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, Qold, Qnew);

            // ---------------------------------------------------------
            // Take a full time step of size dt
            switch( global_ini_params.get_time_order() )
            {


                case 4:

                SetBndValues(Qnew);

                // -- Stage 1 -- //
                ConstructIntegratedR( 0.5*dt, Qnew, smax, F, G);

                // That call is equivalent to the following call:
                // Note that the dt has been rescaled in order to retain the
                // correct units for the flux splitting that will occur in a
                // second.
//              ConstructIntegratedR( 0.5*dt, 
//                  1.0, 0.5, aux,     qnew, 
//                  0.0, 0.0, auxstar, qstar,
//                  smax, F, G);

                // Update the solution:
                ConstructLxWL( Qnew, F, G, Lstar, smax);
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + 0.5*dt*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }
                Qstar.set_t( Qnew.get_t() + 0.5*dt );

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //
                ConstructIntegratedR( dt, 
                    1.0, (1.0/6.0), Qnew, 
                    0.0, (1.0/3.0), Qstar,
                    smax, F, G);
                ConstructLxWL( Qstar, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qold.get_t() + dt );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                break;

                case 5:

                // Coeffients chosen to optimize region of absolute stability 
                // along the imaginary axis.
                //
                // rho = 8.209945182837015e-02 chosen to maximize range of 
                //                             absolute stability region

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 1 -- //
                ConstructIntegratedR( 2.0/5.0*dt, 
                    1.0, 0.5, 125./8.*8.209945182837015e-02, Qnew, 
                    0.0, 0.0, 0.0,                           Qstar,
                    smax, F, G);

                ConstructLxWL( Qnew, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qnew.vget( k ) + (2.0/5.0*dt)*Lstar.vget(k);
                    qstar.vset(k, tmp );
                }
                Qstar.set_t( Qnew.get_t() + (2.0/5.0*dt) );

                // Perform any extra work required:
                AfterStep(dt, Qstar );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // -- Stage 2 -- //
                ConstructIntegratedR( dt, 
                    1.0, 0.5, (1.0/16.0), Qnew, 
                    0.0, 0.0, (5.0/48.0), Qstar,
                    smax, F, G);
                ConstructLxWL( Qstar, F, G, Lstar, smax);

                // Update the solution:
#pragma omp parallel for
                for( int k=0; k < numel; k++ )
                {
                    double tmp = qold.vget( k ) + dt*Lstar.vget(k);
                    qnew.vset(k, tmp );
                }
                Qnew.set_t( Qold.get_t() + dt );

                SetBndValues(Qnew);
                SetBndValues(Qstar);

                // Perform any extra work required:
                AfterStep(dt, Qnew );

                break;

                default:
                printf("Error.  Time order %d not implemented for multiderivative\n", global_ini_params.get_time_order() );
                exit(1);

            }

            // do any extra work
            AfterFullTimeStep(dt, Qold, Qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveMD2D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)       // accept
            { m_accept = 1; }
            else                    //reject
            {   
                t = told;
                if( global_ini_params.get_verbosity() )
                {
                    cout<<"FinSolveMD2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                Qnew.copyfrom( Qold  );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln(Qnew);

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

}
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "RKinfo.h"             // Coefficients for the RK method
#include "FinSolveRK.h"         // Functions directly called from this function
#include "app_defined.h"
#include "StateVars.h"          // Information for state variables
#include "IniParams.h"          // Global parameters accessor
using namespace std;

void FinSolveRK( StateVars& Qnew, double tend, double dtv[] )
{

    dTensorBC3& qnew = Qnew.ref_q  ();
    dTensorBC3&  aux = Qnew.ref_aux();

    // Time stepping information
    const double CFL_max      = global_ini_params.get_max_cfl();      // max CFL number
    const double CFL_target   = global_ini_params.get_desired_cfl();  // target CFL number
    double t                  = Qnew.get_t();
    double dt                 = dtv[1];   // Start with time step from last frame
    double cfl                = 0.0;      // current CFL number
    double dtmin              = dt;       // Counters for max and min time step taken
    double dtmax              = dt;

    // Declare information about the Runge-Kutta method
    const int time_order = global_ini_params.get_time_order();
    RKinfo rk;
    SetRKinfo(time_order, rk);
    double tmp_t = 0.;                    // used for fourth-order stepping

    // Grid information
    const int mx   = global_ini_params.get_mx();
    const int my   = global_ini_params.get_my();
    const int meqn = global_ini_params.get_meqn();
    const int maux = global_ini_params.get_maux();
    const int mbc  = global_ini_params.get_mbc();
    const int numel = qnew.numel();

    // Maximum wave speed at each flux interface value.  Note that the size of
    // this tensor is one longer in each direction.
    dTensorBC3 smax( mx+1, my+1, 2, mbc );           

    // Allocate storage for this solver
    StateVars Qold( t, mx, my, meqn, maux, mbc );
    dTensorBC3& qold   = Qold.ref_q();
    dTensorBC3& auxold = Qold.ref_aux();
    Qold.copyfrom( Qnew );

    // Intermediate stages
    StateVars Qstar( t, mx, my, meqn, maux, mbc );
    dTensorBC3&   qstar = Qstar.ref_q();
    dTensorBC3& auxstar = Qstar.ref_aux();
    Qstar.copyfrom( Qnew );

    dTensorBC3   Lstar(mx, my, meqn, mbc);   // Right hand side of ODE
    dTensorBC3    Lold(mx, my, meqn, mbc);

    // Local storage (for 4th- and 5th-order time stepping)
    StateVars    Q1( t, mx, my, meqn, maux, mbc );
    dTensorBC3&  q1   = Q1.ref_q();
    dTensorBC3&  aux1 = Q1.ref_aux();
    Q1.copyfrom( Qnew );

    StateVars    Q2( t, mx, my, meqn, maux, mbc );
    dTensorBC3&  q2   = Q2.ref_q();
    dTensorBC3&  aux2 = Q2.ref_aux();
    Q2.copyfrom( Qnew );

    // ---------------------------------------------- //
    // -- MAIN TIME STEPPING LOOP (for this frame) -- //
    // ---------------------------------------------- //
    int n_step   = 0;                           // Number of time steps taken
    const int nv = global_ini_params.get_nv();  // Maximum allowable time steps
    while( t<tend )
    {
        // initialize time step
        int m_accept = 0;      
        n_step       = n_step + 1;

        // check if max number of time steps exceeded
        if( n_step > nv )
        {
            printf(" Error in FinSolveRK.cpp: "         );
            printf("Exceeded allowed # of time steps \n");
            printf("    n_step = %d\n", n_step          );
            printf("        nv = %d\n", nv              );
            printf("Terminating program.\n"             );
            exit(1);
        }        

        // copy qnew into qold
        Qold.copyfrom( Qnew );

        // keep trying until we get time step that doesn't violate CFL condition
        while( m_accept==0 )
        {

            // set current time
            Qnew.set_t( t );
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // do any extra work
            BeforeFullTimeStep(dt, Qold, Qnew);

            // Take a full time step of size dt
            switch( time_order )
            {

                case 1:  // First order in time

                    // --------------------------------------------------------
                    // Stage #1 (the only one in this case)
                    rk.mstage = 1;
                    SetBndValues( Qnew  );
                    BeforeStep(dt, Qnew );
                    ConstructL(Qnew, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Qnew);
                    AfterStep(dt, Qnew );
                    // --------------------------------------------------------

                    break;

                case 2:  // Second order in time

                    // ---------------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
                    SetBndValues(  Qnew );
                    BeforeStep(dt, Qnew );
                    ConstructL(Qnew, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Qstar);      
                    AfterStep(dt , Qstar );

                    // ---------------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    SetBndValues(Qstar);
                    BeforeStep(dt, Qstar);
                    ConstructL(Qstar, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qstar, Lstar, Qnew);
                    AfterStep(dt, Qnew );
                    // ---------------------------------------------------------

                    break;

                case 3:  // Third order in time  (low-storage SSP method)

                    // Each update looks like:
                    //
                    // qnew = alpha1 * qstar + alpha2 * qnew + beta * dt * L( qstar )


                    // ---------------------------------------------------------
                    // Stage #1
                    //      alpha1 = 1.0
                    //      alpha2 = 0.0
                    //      beta   = 1.0
                    // ---------------------------------------------------------
                    rk.mstage = 1;
                    SetBndValues(Qnew);
                    BeforeStep(dt,Qnew);    
                    ConstructL(Qnew,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Qstar);
                    AfterStep(dt, Qstar);


                    // ---------------------------------------------------------
                    // Stage #2
                    //      alpha1 = 0.75
                    //      alpha2 = 0.25
                    //      beta   = 0.25
                    // ---------------------------------------------------------

                    rk.mstage = 2;
                    SetBndValues(Qstar);
                    BeforeStep(dt, Qstar);
                    ConstructL(Qstar, Lstar, smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qnew, Lstar, Qstar);
                    AfterStep(dt, Qstar);

                    // ---------------------------------------------------------
                    // Stage #3
                    //      alpha1 = 2/3
                    //      alpha2 = 1/3
                    //      beta   = 2/3
                    // ---------------------------------------------------------

                    rk.mstage = 3;
                    SetBndValues(Qstar);
                    BeforeStep(dt,Qstar);
                    ConstructL(Qstar,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Qstar, Lstar, Qnew);   
                    AfterStep(dt, Qnew);
                    // ---------------------------------------------------------

                    break;

                case 4: // Fourth order in time (10-stages) See Pseudocode 3 in
                        //
                        // "Highly Efficient Strong Stability Preserving Runge-Kutta Methods with
                        // Low-Storage Implementations," David I. Ketcheson, SIAM Journal on Scientific 
                        // Computing, 30(4):2113-2136 (2008)
                        //

                    // -----------------------------------------------
                    Q1.copyfrom( Qnew );
                    Q2.copyfrom( Qnew );

                    // Stage: 1,2,3,4, and 5
                    for (int s=1; s<=5; s++)
                    {
                        rk.mstage = s;
                        SetBndValues(Q1);
                        BeforeStep(dt, Q1);
                        ConstructL(Q1, Lstar, smax);
                        if (s==1)
                        {  Lold.copyfrom( Lstar ); }
                        UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                                rk.beta->get(rk.mstage), dt, Q1, Lstar, Q1);
                        AfterStep(dt, Q1);
                    }

                    // Temporary storage
                    #pragma omp parallel for
                    for( int k=0; k < numel; k++ )
                    {
                        double tmp = (q2.vget(k) + 9.0*q1.vget(k))/25.0;
                        q2.vset(k, tmp );
                        q1.vset(k, 15.0*tmp - 5.0*q1.vget(k) );

                    }

                    // Swap the time values as well (L=1)
                    tmp_t = (Q2.get_t() + 9.0*Q1.get_t())/25.0;
                    Q2.set_t( tmp_t );
                    Q1.set_t( 15.0*tmp_t - 5.0*Q1.get_t() );

                    // Stage: 6,7,8, and 9
                    for (int s=6; s<=9; s++)
                    {
                        rk.mstage = s;
                        SetBndValues(Q1);
                        BeforeStep(dt, Q1);
                        ConstructL(Q1, Lstar, smax);
                        UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                                rk.beta->get(rk.mstage), dt, Q1, Lstar, Q1);
                        AfterStep(dt, Q1);
                    }

                    // Stage: 10
                    rk.mstage = 10;
                    SetBndValues(Q1);
                    BeforeStep(dt,Q1);
                    ConstructL(Q1,Lstar,smax);
                    UpdateSoln(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Q2, Lstar, Q1);
                    AfterStep(dt, Q1);

                    Qnew.copyfrom( Q1 );

                    break;

                case 5: // Fifth order in time (8-stages)
                        // TODO - what paper did these coefficients come from?

//                  Q1.copyfrom( Qnew );   // we can remove two replacements
                                           // here
                    q2.setall(0.);      
                    Q2.set_t( 0.);

                    for (int s=1; s<=8; s++)
                    {
                        rk.mstage = s;
                        SetBndValues( Qnew );
                        BeforeStep(dt, Qnew );
                        ConstructL(Qnew, Lstar, smax);
                        if( s==1 )
                        {  Lold.copyfrom( Lstar ); }

                        UpdateSoln( rk.gamma->get(1,s), rk.gamma->get(2,s), rk.gamma->get(3,s), 
                                rk.delta->get(s), rk.beta->get(s), dt,  Qold, Lstar, Qnew, Q2);

                        AfterStep(dt, Qnew);
                    }

// TODO - the time information for this isn't working correctly.
//printf("Q, t1, t2 = %f, %f, %f \n", qnew.get(1,1), Qnew.get_t(), Q2.get_t() );
//assert_lt( fabs( Qnew.get_t() - t ), 1e-8 );

                    break;

                default:

                    printf("WARNING: torder = %d has not been implemented\n", time_order );
                    break;

            }  // End of switch statement over time-order

            // Do any extra work
            AfterFullTimeStep(dt, Qold, Qnew);

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( global_ini_params.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveRK2D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)       // accept
            { m_accept = 1; }
            else                    //reject
            {   
                t = told;
                if( global_ini_params.get_verbosity() )
                {
                    cout<<"FinSolveRK2D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                Qnew.copyfrom( Qold );
            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln( Qnew );

    } // End of while loop

    // set initial time step for next call to DogSolve:
    dtv[1] = dt;

    DeleteRKInfo(rk);

}
#ifndef _FINSOLVE_RK_
#define _FINSOLVE_RK_

#include <string>
#include "StateVars.h"

// ------------------------------------------------------------
// Function definitions
void ConSoln( const StateVars& Q );

double GetCFL(double dt, double dtmax,
        const dTensorBC3& aux,
        const dTensorBC3& smax);

// These four functions get called in the following order for each stage in
// the Runge-Kutta method:
void BeforeStep(double dt, StateVars& Q );
void ConstructL( StateVars& Q, dTensorBC3& Lstar, dTensorBC3& smax);

// Used for orders 1--4:
void UpdateSoln(double alpha1,double alpha2,double beta,double dt,
    const StateVars& Q1, const dTensorBC3& Lstar, StateVars& Qnew);

// Used for fifth-order stepper:
void UpdateSoln(
    double g1, double g2, double g3, double delta, 
    double beta, double dt,
    const StateVars& Qold, const dTensorBC3& Lstar,
    StateVars& Q1, StateVars& Q2);

void AfterStep(double dt, StateVars& Q );

// Called once before each full time step (c.f. BeforeStep, which is called
// before each stage in an RK method)
void BeforeFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew);
void AfterFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew);

// ------------------------------------------------------------
// Runge-Kutta information
void SetRKinfo(int method2, RKinfo& rk);
void DeleteRKInfo(RKinfo& rk);
// ------------------------------------------------------------

#endif
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "StateVars.h"

using namespace std;

void FinSolveUser( StateVars& Qnew, double tend, double dtv[] )
{


/* TODO - write a template here, and pull it from FinSolveRK.  See the
 * template already created in the 1D code.  (-DS)

    while (t<tend)
    {
        // initialize time step
        m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            cout << " Error in DogSolveUser.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold
        CopyQ(qnew,qold);

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            for (j=1-mbc; j<=(melems+mbc); j++)
            { smax.set(j, 0.0e0 ); }

            // ----------------------------------------------------------------
            //
            //    THIS IS WHERE THE USER-DEFINED TIME-STEPPING SCHEME
            //    SHOULD BE ADDED. IN THE DEFAULT FILE: DogSolveUser.cpp,
            //    THE PROGRAM WILL NOW RETURN AN ERROR MESSAGE.
            // 
            // ----------------------------------------------------------------
            cout << endl;
            cout << " No user-defined time-stepping scheme has been defined yet. " << endl;
            cout << " Copy $FINESS/lib/1d/DogSolveUser.cpp into the current " << endl;
            cout << " directory and modify as needed." << endl << endl;
            exit(1);
            // ----------------------------------------------------------------

            // do any extra work      
            AfterFullTimeStep(dt,node,prim_vol,auxstar,aux,qold,qnew);

            // compute cfl number
            cfl = GetCFL(dt,dtv[2],prim_vol,method,aux,smax);

            // output time step information
            if (method[4]>0) 
            {
                cout << setprecision(3);
                cout << "DogSolve1D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
                // accept
            { m_accept = 1; }
            else 
                //reject
            {   
                t = told;
                if (method[4]>0)
                {
                    cout<<"DogSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                CopyQ(qold,qnew);
            }

        }

        // compute conservation and print to file
        ConSoln(method,node,aux,qnew,t,outputdir);

    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

*/

}
#include <iostream>
#include "dog_math.h"
#include "tensors.h"
#include "stdlib.h"

#include "IniParams.h"
using namespace std;

// Compute the maximum CFL number.  This routine assumes that you have saved
// each wave speed in smax.
//
// Input:
// ------
//
//      dt    - size of current time step
//      dtmax - maximum allowable time step
//      aux   - auxiliary array
//      smax  - maximum wave speed.  smax( 1:mx, 1:my, 1:ndim ).
//
//              smax(i,j,1) = max( f'( q_ij ) ), and 
//              smax(i,j,2) = max( g'( q_ij ) ).  (TODO - this description isn't
//              necessarily correct. smax1 = f'( q_{i-1/2, j } ), and
//                                   smax2 = f'( q_{i, j-1/2 } ).
//
// Returns:
// --------
//
//     cfl - the cfl number
//
// See also: SetMaxWaveSpd (application required) and GlobalWaveSpd.
double GetCFL(double dt, double dtmax,
        const dTensorBC3& aux,
        const dTensorBC3& smax)
{

    const int mx    = global_ini_params.get_mx();
    const int my    = global_ini_params.get_my();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

    double cfl   = -100.0;

    if( dt > dtmax )
    {
        cout << endl;
        cout << " Error: dt is out of bounds ... " << endl;
        cout << "     dt = " << dt << endl;
        cout << endl;
        exit(1);
    }

// This can't be done in parallel!!!  (-DS)
//#pragma omp parallel for
    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
    {
        cfl = Max( dt*smax.get(i, j, 1) / dx,  cfl );
        cfl = Max( dt*smax.get(i, j, 2) / dy,  cfl );
    }

    if( cfl > 1.0e8 )
    {
        cout << endl;
        cout << " Error: CFL number is out of bounds ... " << endl;
        cout << "     CFL = " << cfl << endl;
        cout << endl;
        exit(1);
    }

    return cfl;
}

// Alternative method of pulling the CFL number.  This relies on the fact that
// the "maximum" wave speed has already been computed, and saved in 
// alpha1, and alpha2.
double GetCFL(double dt, double dtmax, double alpha1, double alpha2 )
{

    if( dt > dtmax )
    {
        cout << endl;
        cout << " Error: dt is out of bounds ... " << endl;
        cout << "     dt = " << dt << endl;
        cout << endl;
        exit(1);
    }

    double cfl = -100.0;
    cfl = Max( dt*alpha1 / global_ini_params.get_dx(),  cfl );
    cfl = Max( dt*alpha2 / global_ini_params.get_dy(),  cfl );

    if( cfl > 1.0e8 )
    {
        cout << endl;
        cout << " Error: CFL number is out of bounds ... " << endl;
        cout << "     CFL = " << cfl << endl;
        cout << endl;
        exit(1);
    }

    return cfl;
}
#include "dogdefs.h"
#include "dog_math.h"
#include "GlobalWaveSpd.h"
#include "StateVars.h"
#include "IniParams.h"

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
    const StateVars& Q, 
    double& alpha1, double& alpha2)
{

    const dTensorBC3& q   = Q.const_ref_q  ();
    const dTensorBC3& aux = Q.const_ref_aux();

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    // Needed to define derivatives
    const double dx    = global_ini_params.get_dx();
    const double dy    = global_ini_params.get_dy();
    const double xlow  = global_ini_params.get_xlow();
    const double ylow  = global_ini_params.get_ylow();

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

        double s1, s2;
        SetWaveSpd(nvecx, xedge, Ql, Ql, Auxl, Auxl, s1, s2 );
        alpha1 = Max( alpha1, Max( fabs(s1), fabs(s2) ) );

        SetWaveSpd(nvecy, xedge, Ql, Ql, Auxl, Auxl, s1, s2 );
        alpha2 = Max( alpha2, Max( fabs(s1), fabs(s2) ) );

    }

}
#ifndef _GLOBAL_WAVE_SPD_H_
#define _GLOBAL_WAVE_SPD_H_

// Compute a global wave speed.

// Note: each Numerical flux should be consistent.  That is, if Ql = Qr, then
// fhat( Ql, Qr ) = f( Ql ), where f is the real flux.  For now, we will rely on
// this fact in order to construct the correct maximum wave speeds.
void SetWaveSpd(const dTensor1& nvec, 
        const dTensor1& xedge,
        const dTensor1& Ql, 
        const dTensor1& Qr,
        const dTensor1& Auxl, 
        const dTensor1& Auxr,
        double& s1,double& s2);

#endif
#include "tensors.h"
#include "IniParams.h"

void GridSetup( dTensor3& node, dTensor2& prim_vol)
{

    const double xlow = global_ini_params.get_xlow();
    const double dx   = global_ini_params.get_dx();

    const double ylow = global_ini_params.get_ylow();
    const double dy   = global_ini_params.get_dy();

// TODO - ??? - mnodes = mx*my, or mnodes = (mx,my) ???
//  const int        = node.getsize(1);
//  const int melems = prim_vol.getsize();

//      // Set variable "node"
//  #pragma omp parallel for
//      for (int i=1; i<=mnodes; i++)
//      {  node.set(i,1, xlow+(double(i)-1.0e0)*dx );  }

//      // Set grid cell volumes
//  #pragma omp parallel for
//      for (int j=1; j<=melems; j++)
//      {  prim_vol.set(j, dx);  }

}
#include "tensors.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
//     Dummy Function Call if not using LaxWendroff
//
void D2FluxFunc(const dTensor2& xpts, 
		const dTensor2& Q,
		const dTensor2& Aux,
		dTensor5& D2flux)
{

    const int numpts=xpts.getsize(2);
    D2flux.setall(0.);

}
#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// The expected format is Dflux.get(:,i,j,1:2) = 
//            (\partial f_i, \partial q_j, \partial g_i, \partial q_j )
//
//     Simple advection equation, f'(q) = u1, g'(q) = u2.
//     Burger's equation, f'(q) = q, g'(q) = q
//
void DFluxFunc(const dTensor2& xpts, 
	       const dTensor2& Q,
	       const dTensor2& Aux,
	       dTensor4& Dflux)
{

    const int numpts = xpts.getsize(1);
    for(int i=1; i<=numpts; i++)
    {

        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        Dflux.set(i,1,1,1, 0. );  // First component of flux func
        Dflux.set(i,1,1,2, 0. );  // Second component of flux func
    }

}
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"
using namespace std;

void Output( const StateVars& Q, int nframe )
{

    const dTensorBC3& q   = Q.const_ref_q  ();
    const dTensorBC3& aux = Q.const_ref_aux();
    const double t        = Q.get_t();

    const int meqn    = global_ini_params.get_meqn();
    const int maux    = global_ini_params.get_maux();
    const int mx      = global_ini_params.get_mx();
    const int my      = global_ini_params.get_my();

    string outputdir = global_ini_params.get_output_dir();

    // Open file -- q
    ostringstream fname1;
    fname1 << outputdir << "/" << "q" << setfill('0') 
        << setw(4) << nframe << ".dat";
    ofstream q_file(fname1.str().c_str(), ios::out );

    q_file << setprecision(16);
    q_file << setw(24) << scientific << t << endl;

    // Output each coefficient - TODO, we could potentially reverse the order
    // here, but then the plotting routines will have to change as well.  It
    // is faster to index the arrays in odometer order. (-DS).
    for (int m=1; m<=meqn; m++)
    for (int j=1; j<=my; j++)      
    for (int i=1; i<=mx; i++)      
    {
        q_file << setw(24) << scientific << q.get(i,j,m) << endl;
    }
    q_file.close();

    // Open file -- aux
    ostringstream fname2;
    fname2 << outputdir << "/" << "a" << setfill('0') 
        << setw(4) << nframe << ".dat";
    ofstream aux_file(fname2.str().c_str(), ios::out );

    aux_file << setprecision(16);
    aux_file << setw(24) << scientific << t << endl;

    // Output aux array
    for (int m=1; m<=maux; m++)
    for (int j=1; j<=my; j++)      
    for (int i=1; i<=mx; i++)      
    {
        aux_file << setw(24) << scientific << aux.get(i,j,m) << endl;
    }
    aux_file.close();

    // Output additional information if needed - TODO reintroduce this call
//  void Output_Extra(const dTensor2& node, 
//          const dTensorBC2& aux,
//          const dTensorBC2& q,
//          double t,
//          int nframe,
//          string outputdir);
//  Output_Extra(node,aux,q,t,nframe,outputdir);

}
#ifndef _RKINFO_H_
#define _RKINFO_H_

// Runge-Kutta information
struct RKinfo
{
  int mstage;
  int num_stages;

  dTensor1* alpha1;
  dTensor1* alpha2;
  dTensor1* beta;

  // These two are needed for 5th order Stepping
  dTensor2* gamma;
  dTensor1* delta;

};

#endif
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "dogdefs.h"
#include "IniParams.h"            // accessors for the parameters.ini file
#include "RunFinpack.h"           // Function declarations
#include "StateVars.h"

/*
 * Top level function to run FINESS.
 *
 * See also: $FINESS/lib/main_global.cpp, $FINESS/lib/[1-3]d/RunFinpack.cpp.
 *
 */
int RunFinpack( )
{

    using std::cout;
    using std::endl;
    using std::string;
    using std::scientific;
    using std::setw;
    using std::setprecision;

    // Output title information
    cout << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << "   | FINESS:  The FINite difference ESSentially   |   " << endl;
    cout << "   |          non-oscillatory software Package    |   " << endl;
    cout << "   | Developed by the research group of           |   " << endl;
    cout << "   |            David C. Seal                     |   " << endl;
    cout << "   |            Department of Mathematics         |   " << endl;
    cout << "   |            Michigan State University         |   " << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << endl;


    // Get parameters and print to screen
    cout << global_ini_params.ini_doc_as_string() << endl;
    const IniParams::TimeSteppingMethod::enum_type time_stepping_method = 
	  global_ini_params.get_time_stepping_method();

    // Print information about the parameters to file.  In order to use the
    // MATLAB plotting routines, this call is necessary to pull information
    // from the parameters.ini file.
    WriteQhelp( );

    const int&     nout     = global_ini_params.get_nout();
    const double&  tfinal   = global_ini_params.get_tfinal();
    double dtv[2+1];
    dtv[1] = global_ini_params.get_initial_dt();
    dtv[2] = global_ini_params.get_max_dt();
    const int&     meqn     = global_ini_params.get_meqn();
    const int&     maux     = global_ini_params.get_maux();
    const int&     mx       = global_ini_params.get_mx();
    const int&     my       = global_ini_params.get_my();
    const int&     mbc      = global_ini_params.get_mbc();

    // Dimension arrays
    StateVars Qnew(0., mx, my, meqn, maux, mbc );
    dTensorBC3& qnew = Qnew.ref_q();
    dTensorBC3& aux  = Qnew.ref_aux();

    // Set any auxiliary variables on computational grid
    if( maux > 0 )
    {  SampleFunction(1-mbc, mx+mbc, 1-mbc, my+mbc, qnew, aux, aux, &AuxFunc);  }

    // Set initial data on computational grid
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, qnew, aux, qnew, &QinitFunc);

    // Run AfterQinit to set any necessary variables
    AfterQinit( Qnew );

    // Output initial data to file
    Output( Qnew, 0 );

    // Compute conservation and print to file
    ConSoln( Qnew );

    // Main loop for time stepping
    double tstart = 0.0;
    double tend   = 0.0;
    double dtout  = tfinal/double(nout);    
    for (int n=1; n<=nout; n++)
    {        

        tstart = tend;      assert_lt( fabs(Qnew.get_t()-tend), 1e-13 );
        tend   = tstart + dtout;

        // Solve hyperbolic system from tstart to tend
        if (time_stepping_method == IniParams::TimeSteppingMethod::RK)
        {  
            // Runge-Kutta time-stepping scheme
            FinSolveRK( Qnew, tend, dtv );
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::LxW)
        {
            // Lax-Wendroff time stepping
            FinSolveLxW(Qnew, tend, dtv );
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::MD)
        {
            // Multiderivative time stepping
            FinSolveMD(Qnew, tend, dtv );
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::USER_DEFINED)
        {
            // User-defined time-stepping scheme
            FinSolveUser(Qnew, tend, dtv );
        }
        else
        {
            printf("Time stepping method not implemented\n");
            exit(1);
        }

        // Output data to file
        Output( Qnew, n );

        // Done with solution from tstart to tend
        cout << setprecision(5);
        cout << "FINESS: Frame " << setw(3) << n;
        cout << ": plot files done at time t =";
        cout << setw(12) << scientific << tend << endl;
        cout << endl;
    }

    return 0;
}

// Wrapper functions to make the calls to Qinit and AuxFunc make sense when
// passed into SampleFunction
void QinitFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals)
{
    void QinitFunc(const dTensor2& xpts, dTensor2& qvals);
    QinitFunc(xpts,qvals);
}

void AuxFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals)
{
    void AuxFunc(const dTensor2& xpts, dTensor2& auxvals);
    AuxFunc(xpts, auxvals);
}
#ifndef _RUN_FINESS_H_
#define _RUN_FINESS_H_

#include "StateVars.h"

void WriteQhelp( void );

void Output( const StateVars& Qnew, int nframe );
void QinitFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals);
void AuxFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals);
void AfterQinit( StateVars& Qnew );

void ConSoln( const StateVars& Q );

// Time stepping methods:
void FinSolveRK     ( StateVars& Qnew, double tend, double dtv[] );
void FinSolveLxW    ( StateVars& Qnew, double tend, double dtv[] );
void FinSolveMD     ( StateVars& Qnew, double tend, double dtv[] );
void FinSolveUser   ( StateVars& Qnew, double tend, double dtv[] );

void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC3& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

// ------------------------------------------------------------

#endif
#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "stdlib.h"


#include "IniParams.h"
using namespace std;

// All-purpose routine for sampling a function, and saving its data into a
// single tensor.
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensorBC3 auxin(1-mbc:mx+mbc, 1-mbc:my+mbc, maux    )
//           dTensorBC3   qin(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn    )
//           dTensorBC3  Fout(1-mbc:mx+mbc, 1-mbc:my+mbc, mlength )
// ---------------------------------------------------------------------
//
// The reason there is an extra awkward parameter, mpoints in here is to keep
// the same user interface that DoGPack uses.
//
void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC3& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&))
{    

    const int meqn    = global_ini_params.get_meqn();
    const int maux    = global_ini_params.get_maux();
    const int mlength = Fout.getsize(3);

    // -----------------
    // Quick error check
    // -----------------
    if( meqn<1 || maux <0 || mlength<1 )
    {
        cout << " Error in SampleFunction.cpp ... " << endl;
        cout << "         meqn = " << meqn << endl;
        cout << "         maux = " << maux << endl;
        cout << "      mlength = " << mlength << endl;
        cout << "       istart = " << istart << endl;
        cout << "         iend = " << iend << endl;
        cout << endl;
        exit(1);
    }

    // Grid information
    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();

    // ----------------------------------
    // Loop over all elements of interest
    // ----------------------------------    

    const int mpoints = 1;
#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    {

        const double xc = xlow + (double(i)-0.5)*dx;
        for (int j=jstart; j<=jend; j++)
        {
            const double yc = ylow + (double(j)-0.5)*dy;

            // each of these three items needs to be private to each thread ..
            dTensor2 xpts(mpoints,2);
            dTensor2 qvals(mpoints, meqn), auxvals(mpoints, maux);
            dTensor2 fvals(mpoints, mlength);

            qvals.setall(0.);
            auxvals.setall(0.);

            // Loop over each quadrature point
            for (int m=1; m<= mpoints; m++)
            {

                // grid point x
                xpts.set( m, 1, xc );
                xpts.set( m, 2, yc );

                // Solution values (q) at each grid point
                for (int me=1; me<=meqn; me++)
                {

                    for (int k=1; k<=mpoints; k++)
                    {
                        qvals.set(m, me, qvals.get(m,me) + qin.get(i, j, me) );
                    }
                }

                // Auxiliary values (aux) at each grid point
                for (int ma=1; ma<=maux; ma++)
                {
                    auxvals.set(m,ma, 0.0 );
                    for (int k=1; k<=mpoints; k++)
                    {
                        auxvals.set(m,ma, auxvals.get(m,ma) + auxin.get(i,j,ma) );
                    }
                }
            }

            // Call user-supplied function to set fvals
            Func(xpts, qvals, auxvals, fvals);

            // Evaluate "integrals"
            for (int m1=1; m1<=mlength; m1++)
            {
                Fout.set(i, j, m1, fvals.get(1,  m1) );
            }

        }
    }

}

// All-purpose routine for sampling a function, and saving its data into a
// single tensor.
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensorBC3 auxin(1-mbc:mx+mbc, 1-mbc:my+mbc, maux          )
//           dTensorBC3   qin(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn          )
//           dTensorBC4  Fout(1-mbc:mx+mbc, 1-mbc:my+mbc, mlength, ndim )
// ---------------------------------------------------------------------
//
// The reason there is an extra awkward parameter, mpoints in here is to keep
// the same user interface that DoGPack uses.
//
void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC4& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&))
{    

    const int meqn    = global_ini_params.get_meqn();
    const int maux    = global_ini_params.get_maux();
    const int mlength = Fout.getsize(3);

    // -----------------
    // Quick error check
    // -----------------
    if( meqn<1 || maux <0 || mlength<1 )
    {
        cout << " Error in SampleFunction.cpp ... " << endl;
        cout << "         meqn = " << meqn << endl;
        cout << "         maux = " << maux << endl;
        cout << "      mlength = " << mlength << endl;
        cout << "       istart = " << istart << endl;
        cout << "         iend = " << iend << endl;
        cout << endl;
        exit(1);
    }

    // Grid information
    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();

    // ----------------------------------
    // Loop over all elements of interest
    // ----------------------------------    

    const int mpoints = 1;
    const int ndim    = 2;  // length of the vector in Fout
#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    {

        const double xc = xlow + (double(i)-0.5)*dx;
        for (int j=jstart; j<=jend; j++)
        {
            const double yc = ylow + (double(j)-0.5)*dy;

            // each of these three items needs to be private to each thread ..
            dTensor2 xpts(mpoints,2);
            dTensor2 qvals(mpoints, meqn), auxvals(mpoints, maux);
            dTensor3 fvals(mpoints, mlength, ndim);

            qvals.setall(0.);
            auxvals.setall(0.);

            // Loop over each quadrature point
            for (int m=1; m<= mpoints; m++)
            {

                // grid point x
                xpts.set( m, 1, xc );
                xpts.set( m, 2, yc );

                // Solution values (q) at each grid point
                for (int me=1; me<=meqn; me++)
                {

                    for (int k=1; k<=mpoints; k++)
                    {
                        qvals.set(m, me, qvals.get(m,me) + qin.get(i, j, me) );
                    }
                }

                // Auxiliary values (aux) at each grid point
                for (int ma=1; ma<=maux; ma++)
                {
                    auxvals.set(m,ma, 0.0 );
                    for (int k=1; k<=mpoints; k++)
                    {
                        auxvals.set(m,ma, auxvals.get(m,ma) + auxin.get(i,j,ma) );
                    }
                }
            }

            // Call user-supplied function to set fvals
            Func(xpts, qvals, auxvals, fvals);

            // Evaluate "integrals"
            for(int m1=1; m1<=mlength; m1++)
            for(int  d=1;  d<=ndim;     d++)
            {
                Fout.set(i, j, m1, d, fvals.get(1,m1,d) );
            }

        }
    }

}
#include "tensors.h"
#include "RKinfo.h"

void SetRKinfo(int method2, RKinfo& rk)
{
    rk.mstage = 0;

    switch( method2 )
    {
        case 5:
            rk.num_stages = 8;
            break;
        case 4:
            rk.num_stages = 10;
            break;
        default:
            rk.num_stages = method2;
    }

    rk.alpha1 = new dTensor1(rk.num_stages);
    rk.alpha2 = new dTensor1(rk.num_stages);
    rk.beta   = new dTensor1(rk.num_stages);

    // 5th order stuff 
    rk.delta  = new dTensor1(rk.num_stages);
    rk.gamma  = new dTensor2(3, rk.num_stages);

    switch(method2)
    {
        case 1: // first-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            break;

        case 2: // second-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            rk.alpha1->set(2, 0.5 );
            rk.alpha2->set(2, 0.5 );
            rk.beta->set(  2, 0.5 );

            break;

        case 3: // third-order

            rk.alpha1->set(1, 1.0 );
            rk.alpha2->set(1, 0.0 );
            rk.beta->set(  1, 1.0 );

            rk.alpha1->set(2, 0.75 );
            rk.alpha2->set(2, 0.25 );
            rk.beta->set(  2, 0.25 );

            rk.alpha1->set(3, 2.0e0/3.0e0 );
            rk.alpha2->set(3, 1.0e0/3.0e0 );
            rk.beta->set(  3, 2.0e0/3.0e0 );

            break;

        case 4: // fourth-order

            for (int i=1; i<=9; i++)
            {
                rk.alpha1->set(i, 1.0 );
                rk.alpha2->set(i, 0.0 );
                rk.beta->set(  i, 1.0/6.0 );
            }

            rk.alpha1->set(10, 1.0 );
            rk.alpha2->set(10, 3.0/5.0 );
            rk.beta->set(  10, 1.0/10.0 );

            break;

        case 5: // 5th-order

            rk.delta->set(1, 1.0e0);
            rk.delta->set(2, 1.528486658778845e00);
            rk.delta->set(3, 4.720094096662784e-02);
            rk.delta->set(4, 8.801244253465348e-01);
            rk.delta->set(5, 1.019066090228480e+00);
            rk.delta->set(6, 1.049772291176110e+01);
            rk.delta->set(7, -4.254616508506826e+00);
            rk.delta->set(8, 0.0);

            rk.gamma->set(1, 1, 0.0);
            rk.gamma->set(1, 2, -1.552288007713033e+01);
            rk.gamma->set(1, 3,  4.127286635722417e-01);
            rk.gamma->set(1, 4, -1.011819196331377e+00);
            rk.gamma->set(1, 5, -2.765748383780848e-01);
            rk.gamma->set(1, 6,  5.075770311217778e-02);
            rk.gamma->set(1, 7,  6.999810478513669e+00);
            rk.gamma->set(1, 8, -1.114908881433104e+01);

            rk.gamma->set(2, 1,  1.0);
            rk.gamma->set(2, 2,  6.534691420958578e+00);
            rk.gamma->set(2, 3,  2.280056542904473e-01);
            rk.gamma->set(2, 4,  1.308684311397668e+00);
            rk.gamma->set(2, 5,  4.769419552531064e-01);
            rk.gamma->set(2, 6, -6.368809762042849e-03);
            rk.gamma->set(2, 7,  9.339446057238532e-02);
            rk.gamma->set(2, 8,  9.556626047962331e-01);

            rk.gamma->set(3, 1, 0.0);
            rk.gamma->set(3, 2, 0.0);
            rk.gamma->set(3, 3, 0.0);
            rk.gamma->set(3, 4, -2.510747784045939e+00);
            rk.gamma->set(3, 5, -8.576822794622042e-01);
            rk.gamma->set(3, 6,  1.044599944472272e+00);
            rk.gamma->set(3, 7, -7.000810861049136e+00);
            rk.gamma->set(3, 8,  1.906311811144179e+00);

            rk.beta->set(1,  8.653258038183180e-02);
            rk.beta->set(2,  9.544677980851571e-01);
            rk.beta->set(3,  2.651941386774408e-01);
            rk.beta->set(4,  2.736914413910379e-01);
            rk.beta->set(5,  5.999778649323600e-01);
            rk.beta->set(6,  4.433177471748104e-03);
            rk.beta->set(7,  5.309971130968292e-03);
            rk.beta->set(8,  5.830861806762871e-01);
            break;

    }
}

void DeleteRKInfo(RKinfo& rk)
{
    delete rk.alpha1;
    delete rk.alpha2;
    delete rk.beta;
    delete rk.gamma;
    delete rk.delta;
}
#ifndef _STATEVARS_H_
#define _STATEVARS_H_

#include <algorithm>

#include "tensors.h"


class StateVars{
    private:
        double t;
        dTensorBC3 q;
        dTensorBC3 aux;
    public:
        StateVars(double t, int mx, int my, int meqn, int maux, int mbc ):            
            t(t),
            q(mx, my, meqn, mbc),
            aux(mx, my, std::max(1, maux), mbc)
        { }

        double get_t() const{
            return this->t;
        }

        void set_t(double t){
            this->t = t;
        }
        
        const dTensorBC3& const_ref_q() const{
            return this->q;
        }

        dTensorBC3& ref_q(){
            return this->q;
        }

        const dTensorBC3& const_ref_aux() const{
            return this->aux;
        }

        dTensorBC3& ref_aux(){
            return this->aux;
        }

        // Copy the contents from another state variable to this state
        // variable
        void copyfrom( const StateVars& Qin )
        {
            this->q.copyfrom( Qin.const_ref_q() );
            this->aux.copyfrom( Qin.const_ref_aux() );
            this->t = Qin.get_t();
        }
};

#endif
#include "tensors.h"
#include "StateVars.h"

// Update the solution using the constructed Lstar
//
// The Low-Storage RK methods use a combination of two different values of q,
// together with a right hand side, L.  This routine computes the following:
//
//     qnew = alpha1*qstar + alpha2*qnew + beta*dt*Lstar.
//
// The main loop covers all elements (including boundary cells) of the arrays.
//
// See also: RKinfo.h and SetRKinfo.cpp
void UpdateSoln(double alpha1, double alpha2, double beta, double dt,
    const StateVars& Qstar, const dTensorBC3& Lstar, StateVars& Qnew)
{

    const dTensorBC3& qstar   = Qstar.const_ref_q();
    const dTensorBC3& auxstar = Qstar.const_ref_aux();
    double tstar              = Qstar.get_t  ();

    dTensorBC3&  qnew   = Qnew.ref_q  ();
    dTensorBC3& auxnew  = Qnew.ref_aux();
    double tnew         = Qnew.get_t  ();

    // Update time
    Qnew.set_t( alpha1*tstar + alpha2*tnew + beta*dt );

    const int numel = qnew.numel();
#pragma omp parallel for
    for( int k=0; k < numel; k++ )
    {
        double tmp = alpha1*qstar.vget(k)+alpha2*qnew.vget(k)+beta*dt*Lstar.vget(k);
        qnew.vset( k, tmp );
    }


}

// Update the solution using the constructed Lstar
//
// This version of UpdateSoln is used for the fifth-order time stepping.
void UpdateSoln(
    double g1, double g2, double g3, double delta, 
    double beta, double dt,
    const StateVars& Qold, const dTensorBC3& Lstar,
    StateVars& Q1, StateVars& Q2)
{

    const dTensorBC3& qold    = Qold.const_ref_q();
    const dTensorBC3& auxold  = Qold.const_ref_aux();

    dTensorBC3&  q1   = Q1.ref_q  ();
    dTensorBC3& aux1  = Q1.ref_aux();

    dTensorBC3&  q2   = Q2.ref_q  ();
    dTensorBC3& aux2  = Q2.ref_aux();

    const int     mx = q1.getsize(1);
    const int     my = q1.getsize(2);
    const int   meqn = q1.getsize(3);
    const int   maux = aux1.getsize(3);
    const int    mbc = q1.getmbc();

    const int numel = q1.numel();
#pragma omp parallel for
    for( int k=0; k < numel; k++ )
    {

        // update q2
        const double s1 = q1.vget( k );
        const double s2 = q2.vget( k ) + delta*s1;
        q2.vset(k, s2 );

        // update q
        const double s3  = qold.vget( k );
        const double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.vget(k);
        q1.vset(k, tmp );

    }

}
#include <sstream>
#include <fstream>
#include <string>
#include "stdlib.h"

#include "IniParams.h"
#include "StateVars.h"

using namespace std;

// This routine is called once per simulation.  The purpose of this routine is
// to provide access to the plotting routines provided from DoGPack.
//
// Make sure that the global IniParams has already been called before calling
// this function.
//
// See: $FINESS/viz/matlab/plotfin1.m for where this gets parsed in the
// plotting routines.
void WriteQhelp( )
{

    string outputdir = global_ini_params.get_output_dir();

    ostringstream fname;
    fname << outputdir << "/" << "qhelp.dat";

    // Open file -- qhelp.dat
    ostringstream fname1;
    fname1 << outputdir << "/" << "qhelp.dat";
    FILE* file = fopen( fname1.str().c_str(), "w" );

    fprintf(file,"%16d   : ndims   \n", global_ini_params.get_ndims() );
    fprintf(file,"%16d   : meqn    \n", global_ini_params.get_meqn()  );
    fprintf(file,"%16d   : maux    \n", global_ini_params.get_maux()  );
    fprintf(file,"%16d   : nout    \n", global_ini_params.get_nout()  );
    fprintf(file,"%16d   : mx      \n", global_ini_params.get_mx()    );
    fprintf(file,"%16d   : my      \n", global_ini_params.get_my()    );
    fprintf(file,"%16.8e : xlow    \n", global_ini_params.get_xlow()  );
    fprintf(file,"%16.8e : xhigh   \n", global_ini_params.get_xhigh() );
    fprintf(file,"%16.8e : ylow    \n", global_ini_params.get_ylow()  );
    fprintf(file,"%16.8e : yhigh   \n", global_ini_params.get_yhigh() );
    fclose(file);

}
#ifndef _APP_SPECIFIC_H_
#define _APP_SPECIFIC_H_

#include "tensors.h"
#include "StateVars.h"

void QinitFunc(const dTensor2& xpts, dTensor2& qvals);
void AuxFunc(const dTensor2& xpts, dTensor2& auxvals);
void FluxFunc(const dTensor2& xpts, const dTensor2& Q, const dTensor2& Aux, dTensor3& flux);
void SetBndValues( StateVars& Q );
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, const dTensor2& auxvals, dTensor2& source);

#endif
#ifndef _CONSTRUCTS_H_
#define _CONSTRUCTS_H_

#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "assert.h"

void ConstructIntegratedR( double dt, const StateVars& Q,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

void LocalIntegrate( 
    int nterms, double dx, double dy, double xc, double yc,
    int meqn, int maux, int mpts_sten, int half_mpts_sten,
    const int i, const int j, const dTensorBC3& q, const dTensorBC3& aux, 
    const dTensorBC4& R, 
    dTensor1& f_t, dTensor1& f_tt,
    dTensor1& g_t, dTensor1& g_tt
    );

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1,
    const StateVars& Q1,
    double alpha2, double beta2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

void ConstructIntegratedR( double dt, 
    double alpha1, double beta1, double charlie1, 
    const StateVars& Q1,
    double alpha2, double beta2, double charlie2,
    const StateVars& Q2,
    dTensorBC3& smax, dTensorBC3& F, dTensorBC3& G);

void ConstructLxWL( const StateVars& Q,
        dTensorBC3& F,         // <--- new term: integrated flux, f
        dTensorBC3& G,         // <--- new term: integrated flux, g
        dTensorBC3& Lstar,
        dTensorBC3& smax);

#endif
#define NDIMS 2
#ifndef _HOOKS_H_
#define _HOOKS_H_
#include "tensors.h"
#include "StateVars.h"

void AfterQinit( StateVars& Q );
void BeforeFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew );
void BeforeStep(double dt, StateVars& Q);
void AfterStep(double dt, StateVars& Q );
void AfterFullTimeStep(double dt, StateVars& Qold, StateVars& Qnew );

#endif
#ifndef _MISC2D_H_
#define _MISC2D_H_

#include "tensors.h"
#include "StateVars.h"

double GetCFL(double dt, double dtmax, const dTensorBC3& aux, const dTensorBC3& smax);
double GetCFL(double dt, double dtmax, double alpha1, double alpha2 );
void ConSoln( const StateVars& Q );

#endif
#include <string>
#include <stdlib.h> // for system()
#include "assert.h"
#include "debug.h"

#include "IniParams.h"

// Function used to call startscript from $(FINESS)/scripts;
//
// That shell script creates the output directory and copies a few files in
// the output folder
void RunStartScript(std::string parameters_ini_filename)
{
    char command_str[1024];
//    const char* get_outputdir();
    std::string output_dir = global_ini_params.get_output_dir();

    // run startscript
    // to create output directory if it does not exist
    // and copy data files to output directory
    int numchars = snprintf(command_str,1024,
            "if test -f startscript && test -x startscript;\n"
            "then ./startscript %s %s\n"
            "else ${FINESS}/scripts/startscript %s %s\n"
            "fi", output_dir.c_str(),
            parameters_ini_filename.c_str(),
            output_dir.c_str(),
            parameters_ini_filename.c_str());
    assert_lt(numchars,1023);
    assert_gt(numchars,0);
    int exit_status = system(command_str);
}

#include <cmath>
#include <stdexcept>
#include "assert.h"
#include "tensors.h"
#include "WenoReconstruct.h"
#include "IniParams.h"

// All purpose routine for computing a conservative WENO reconstruction that
// is based on matching polynomials with cell averages.  Upon taking
// differences of these values, you get a high-order approximation to the
// derivative.
//
// For example, g_x( x_i ) \approx ( g_{i+1/2} - g_{i-1/2} ) / dx.
//
// Input:
//
//      g( 1:meqn, 1:ws ) - list of meqn functions to be differentiated. ws =
//                          size of stencil under consideration.  ws =
//                          space_order.
//
// Output:
//
//      g_reconst( 1:meqn, 1 ) - The reconstructed value of g evaluated at 
//                          the 'right' half of the stencil, i+1/2.  
//
//      To get the other value at i-1/2, reverse the stencil, and call 
//      this same function again.
//
//
// For example, in WENO5, one passes in the following stencil:
//
//     u = { u_{i-2}, u_{i-1}, u_i, u_{i+1}, u_{i+2} },
//
// and then reconstructs the value u_{i+1/2} with this method.  To get the
// value u_{i-1/2}, pass in the following stencil:
//
//     u_{i-1/2} = { u_{i+2}, u_{i+1}, u_i, u_{i-1}, u_{i-2} }.
//
//void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst)
reconstruct_t GetWenoReconstruct()
{

    // TODO - implement WENO3 and WENO-Z version
    if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 5)
        return &WenoReconstruct_JS5;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 7)
        return &WenoReconstruct_JS7;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 9)
        return &WenoReconstruct_JS9;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 11)
        return &WenoReconstruct_JS11;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::FD && global_ini_params.get_space_order() == 5)
        return &WenoReconstruct_FD5;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::FD && global_ini_params.get_space_order() == 7)
        return &WenoReconstruct_FD7;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::FD && global_ini_params.get_space_order() == 9)
        return &WenoReconstruct_FD9;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::Z  && global_ini_params.get_space_order() == 5)
        return &WenoReconstruct_Z5;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::Z  && global_ini_params.get_space_order() == 7)
        return &WenoReconstruct_Z7;
    else if(global_ini_params.get_space_order() == 1)
        return &WenoReconstructLLF;
    else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::Z  && global_ini_params.get_space_order() == 9)
    {
        printf("Warning: we're not sure these are the correct coefficients for WENOZ-9!\n");
        return &WenoReconstruct_Z9;
    }
    else
        throw(std::logic_error("Requested WENO Reconstruction not implemented."));
}

// ---------------------------------------------------- // 
// SECTION: WENO-JS reconstructions                     //
// ---------------------------------------------------- // 

// Fifth-order Jiang and Shu WENO reconstruction
//
//  See: G.-S. Jiang and C.-W. Shu, "Efficient Implementation of Weighted ENO 
//       Schemes". J. Comput. Phys. 126 (1996), pp. 202-228. 
//
//  or Section 2.2 of the recent SIAM review paper:
//
//       C.-W. Shu, "High Order Weighted Essentially Nonoscillatory Schemes for
//                   Convection Dominated Problems", SIAM Review, Vol. 51, 
//                   No. 1, pp 82--126, (2009).
//
// See also: WenoReconstruct.  This routine requires a stencil of lenght 5 to
// be passed in.
void WenoReconstruct_JS5( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq( g.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    double g0, g1, g2;           
    g0 = 0.1; g1 = 0.6; g2 = 0.3;


    const double eps         = global_ini_params.get_epsilon(); 
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    const int meqn = g.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // Compute smoothness indicators (identical for left/right values):
        beta0 =(13./12.)*pow(uim2-2*uim1+ui,2)+0.25*pow(uim2-4*uim1+3*ui,2);
        beta1 =(13./12.)*pow(uim1-2*ui+uip1,2)+0.25*pow(uim1-uip1,2);
        beta2 =(13./12.)*pow(ui-2*uip1+uip2,2)+0.25*pow(3*ui-4*uip1+uip2,2);
        
        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
        u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
        u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;
        
        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*pow(eps+beta0,-power_param);
        omt1 = g1*pow(eps+beta1,-power_param);
        omt2 = g2*pow(eps+beta2,-power_param);
        omts = omt0+omt1+omt2;

        // # Return 5th-order conservative reconstruction
        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3)/omts );

    }

}

// 7th-order WENO reconstruction
//
// See: D.S. Balsara and C.-W. Shu, "Monotonicity Preserving Weighted 
//      Essentially Non-oscillatory Schemes with Increasingly High Order of 
//      Accuracy". J. Comput. Phys. 160 (2000), pp. 405-452.
//
// See also: WenoReconstruct.  This routine requires a stencil of length 7 to
// be passed in.
void WenoReconstruct_JS7( const dTensor2& g, dTensor2& g_reconst )
{
    assert_eq(g.getsize(2), 7);
    double uim3, uim2, uim1, ui, uip1, uip2, uip3;
    double u1, u2, u3, u4;
    
    double beta0, beta1, beta2, beta3;
    double omt0, omt1, omt2, omt3, omts;

    const double g0 = 1./35., g1 = 12./35., g2 = 18./35., g3 = 4./35.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++){
        uim3 = g.get(m, 1);
        uim2 = g.get(m, 2);
        uim1 = g.get(m, 3);
        ui   = g.get(m, 4);
        uip1 = g.get(m, 5);
        uip2 = g.get(m, 6);
        uip3 = g.get(m, 7);

        beta0 = uim3*(  547.*uim3 -  3882.*uim2  + 4642.*uim1 - 1854.*ui) + 
                uim2*( 7043.*uim2 - 17246.*uim1  + 7042.*ui) +              
                uim1*(11003.*uim1 -  9402.*ui  ) + 2107. * pow(ui, 2);
        beta1 = uim2*(  267.*uim2 - 1642.*uim1   + 1602.*ui - 494.*uip1) + 
                uim1*( 2843.*uim1 - 5966.*ui     + 1922.*uip1) +
                ui*(   3443.*ui   - 2522.*uip1 ) + 547. * pow(uip1, 2);
        beta2 = uim1*( 547.*uim1 - 2522.*ui     + 1922.*uip1 - 494.*uip2) + 
                ui  *(3443.*ui   - 5966.*uip1   + 1602.*uip2 )   + 
                uip1*(2843.*uip1 - 1642.*uip2 ) + 267.*pow(uip2, 2);
        beta3 = ui*  ( 2107.*ui   -  9402.*uip1   + 7042.*uip2 - 1854.*uip3 ) + 
                uip1*(11003.*uip1 - 17246.*uip2   + 4642.*uip3 )              + 
                uip2*( 7043.*uip2 -  3882.*uip3 ) + 547.*pow(uip3, 2);
        u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui;
        u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1;
        u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2;
        u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3;

        omt0 = g0*pow(eps+beta0, -power_param);
        omt1 = g1*pow(eps+beta1, -power_param);
        omt2 = g2*pow(eps+beta2, -power_param);
        omt3 = g3*pow(eps+beta3, -power_param);
        omts = omt0 + omt1 + omt2 + omt3;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4)/omts);
    }
}

// 9th-order WENO reconstruction
//
// See also: WenoReconstruct.  This routine requires a stencil of length 9 to
// be passed in.
void WenoReconstruct_JS9( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 9);
    double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4;
    double u1, u2, u3, u4, u5;
    
    double beta0, beta1, beta2, beta3, beta4;
    double omt0, omt1, omt2, omt3, omt4, omts;

    const double g0 = 1.0/126.0, g1 = 10./63., g2 = 10./21., g3 = 20./63., g4 = 5./126.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++)
    {

        uim4 = g.get(m, 1);
        uim3 = g.get(m, 2);
        uim2 = g.get(m, 3);
        uim1 = g.get(m, 4);
        ui   = g.get(m, 5);
        uip1 = g.get(m, 6);
        uip2 = g.get(m, 7);
        uip3 = g.get(m, 8);
        uip4 = g.get(m, 9);

        beta0 = uim4*(22658.*uim4 - 208501.*uim3 + 364863.*uim2 - 288007.*uim1 + 86329.*ui) +
                uim3*(482963.*uim3 - 1704396.*uim2 + 1358458.*uim1 - 411487. * ui) +
                uim2*(1521393.*uim2 - 2462076.*uim1 + 758823.*ui) + 
                uim1*(1020563.*uim1 - 649501.*ui) +
                107918.*pow(ui, 2);
        beta1 = uim3*(6908.*uim3  - 60871.*uim2 + 99213.*uim1 - 70237.*ui + 18079.*uip1) +
                uim2*(138563.*uim2 - 464976.*uim1 + 337018.*ui - 88297.*uip1) +
                uim1*(406293.*uim1 - 611976*ui + 165153.*uip1) +
                ui*(242723.*ui - 140251.*uip1) + 
                22658.*pow(uip1, 2);
        beta2 = uim2*(6908.*uim2 - 51001.*uim1 + 67923.*ui - 38947.*uip1 + 8209.*uip2) +
                uim1*(104963.*uim1 - 299076.*ui + 179098.*uip1 - 38947.*uip2) +
                ui*(231153.*ui - 299076.*uip1 + 67923.*uip2) +
                uip1*(104963.*uip1 - 51001.*uip2) +
                6908.*pow(uip2, 2);
        beta3 = uim1*(22658.*uim1 - 140251.*ui + 165153.*uip1 - 88297.*uip2 + 18079.*uip3) +
                ui*(242723.*ui - 611976.*uip1 + 337018.*uip2 - 70237.*uip3) +
                uip1*(406293.*uip1 - 464976.*uip2 + 99213.*uip3) +
                uip2*(138563.*uip2 - 60871.*uip3) +
                6908.*pow(uip3, 2);
        beta4 = ui*(107918.*ui - 649501.*uip1 + 758823.*uip2 - 411487.*uip3 + 86329.*uip4) +
                uip1*(1020563.*uip1 - 2462076.*uip2 + 1358458.*uip3 - 288007*uip4) + 
                uip2*(1521393.*uip2 - 1704396.*uip3 + 364863.*uip4) +
                uip3*(482963.*uip3 - 208501.*uip4) +
                22658.*pow(uip4, 2);


        u1 = (1./5.)*uim4   + (-21./20.)*uim3 + (137./60.)*uim2 + (-163./60.)*uim1 + (137./60.)*ui;
        u2 = (-1./20.)*uim3 + (17./60.)*uim2  + (-43./60.)*uim1 + (77./60.)*ui     + (1./5.)*uip1;
        u3 = (1./30.)*uim2  + (-13./60.)*uim1 + (47./60.)*ui    + (9./20.)*uip1    + (-1./20.)*uip2;
        u4 = (-1./20.)*uim1 + (9./20.)*ui     + (47./60.)*uip1  + (-13./60.)*uip2  + (1./30.)*uip3;
        u5 = (1./5.)*ui     + (77./60.)*uip1  + (-43./60.)*uip2 + (17./60.)*uip3   + (-1./20.)*uip4;


        omt0 = g0*pow(eps+beta0, -power_param);
        omt1 = g1*pow(eps+beta1, -power_param);
        omt2 = g2*pow(eps+beta2, -power_param);
        omt3 = g3*pow(eps+beta3, -power_param);
        omt4 = g4*pow(eps+beta4, -power_param);
        omts = omt0 + omt1 + omt2 + omt3 + omt4;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4 + omt4*u5)/omts);
    }
}

// 11th-order WENO reconstruction (Jiang and Shu weights)
void WenoReconstruct_JS11( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 11);
    double uim5, uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4, uip5;
    double u1, u2, u3, u4, u5, u6;
    
    double beta0, beta1, beta2, beta3, beta4, beta5, beta6;
    double omt0, omt1, omt2, omt3, omt4, omt5, omt6, omts;

    const double g0=1.0/462., g1=5./77, g2=25./77., g3=100./231., g4=25./154., g5=1./77.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++)
    {

        uim5 = g.get(m, 1);
        uim4 = g.get(m, 2);
        uim3 = g.get(m, 3);
        uim2 = g.get(m, 4);
        uim1 = g.get(m, 5);
        ui   = g.get(m, 6);
        uip1 = g.get(m, 7);
        uip2 = g.get(m, 8);
        uip3 = g.get(m, 9);
        uip4 = g.get(m,10);
        uip5 = g.get(m,11);

        beta0 = uim5*(  1152561.*uim5 - 12950184. *uim4 + 29442256. *uim3 - 33918804. * uim2 + 19834350.*uim1 - 4712740.*ui ) 
              + uim4*( 36480687.*uim4 - 166461044.*uim3 + 192596472.*uim2 - 113206788.* uim1 + 27060170.*ui )
              + uim3*(190757572.*uim3 - 444003904.*uim2 + 262901672.*uim1 - 63394124. * ui )
              + uim2*(260445372.*uim2 - 311771244.*uim1 + 76206736. *ui )
              + uim1*(94851237. *uim1 - 47460464.*ui )  + 6150211.*ui*ui;

        beta1 = uim4*(  271779.* uim4 - 3015728.*uim3  + 6694608.  * uim2 - 7408908.  *  uim1 + 4067018.*ui - 880548.*uip1 )
              + uim3*( 8449957.* uim3 - 37913324.*uim2 + 42405032. * uim1 - 23510468. * ui    + 5134574.*uip1 )
              + uim2*(43093692.* uim2 - 97838784.*uim1 + 55053752. * ui   - 12183636. * uip1 )
              + uim1*(56662212.* uim1 - 65224244.*ui   + 14742480. * uip1 )
              + ui  *(19365967.* ui   - 9117992.*uip1 ) + 1152561. *uip1*uip1;

        beta2  = uim3*(139633.  *uim3 - 1429976.  *uim2 + 2863984. *uim1 - 2792660.*ui   + 1325006.*uip1 - 245620.*uip2 ) 
               + uim2*(3824847. *uim2 - 15880404. *uim1 + 15929912.*ui   - 7727988.*uip1 + 1458762.*uip2 ) 
               + uim1*(17195652.*uim1 - 35817664. *ui   + 17905032.*uip1 - 3462252.*uip2 ) 
               + ui  *(19510972.*ui   - 20427884. *uip1 + 4086352. *uip2 ) 
               + uip1*(5653317. *uip1 - 2380800.  *uip2 ) + 271779.*uip2*uip2;

        beta3 = 
              uim2*(  271779.* uim2 - 2380800.  *uim1 + 4086352. *ui   - 3462252.*uip1 + 1458762.*uip2 - 245620.* uip3 ) 
            + uim1*(5653317. * uim1 - 20427884. *ui    + 17905032.*uip1 - 7727988.*uip2 + 1325006.*uip3 ) 
            + ui  *(19510972.* ui   - 35817664. *uip1  + 15929912.*uip2 - 2792660.*uip3 ) 
            + uip1*(17195652.* uip1 - 15880404. *uip2  + 2863984. *uip3 ) 
            + uip2*(3824847. * uip2 - 1429976.  *uip3 ) 
            + 139633.*uip3*uip3;

        beta4 = 
            uim1*(1152561. * uim1 - 9117992. * ui   + 14742480.* uip1 - 12183636.*uip2 + 5134574.*uip3 - 880548.*uip4 ) 
          + ui  *(19365967.* ui   - 65224244.* uip1 + 55053752.* uip2 - 23510468.*uip3 + 4067018.*uip4 ) 
          + uip1*(56662212.* uip1 - 97838784.* uip2 + 42405032.* uip3 - 7408908. *uip4 ) 
          + uip2*(43093692.* uip2 - 37913324.* uip3 + 6694608. * uip4 ) 
          + uip3*(8449957.*  uip3 - 3015728. * uip4 ) 
          + 271779.*uip4*uip4;

        beta5 = 
          ui  *(6150211.*   ui   - 47460464. *uip1 + 76206736. *uip2 - 63394124. *uip3 + 27060170.*uip4 - 4712740.*uip5 ) 
        + uip1*(94851237.*  uip1 - 311771244.*uip2 + 262901672.*uip3 - 113206788.*uip4 + 19834350.*uip5 ) 
        + uip2*(260445372.* uip2 - 444003904.*uip3 + 192596472.*uip4 - 33918804. *uip5 ) 
        + uip3*(190757572.* uip3 - 166461044.*uip4 + 29442256. *uip5 ) 
        + uip4*(36480687.*  uip4 - 12950184. *uip5 ) 
        + 1152561.*uip5*uip5;

        u1 = (-1./6.)*uim5 + (31./30.)*uim4 + (-163./60.)*uim3 + (79./20.)*uim2 + (-71./20.)*uim1 + (49./20.)*ui;  
        u2 = (1./30.)*uim4 + (-13./60.)*uim3 + (37./60.)*uim2 + (-21./20.)*uim1 + (29./20.)*ui + (1./6.)*uip1;  
        u3 = (-1./60.)*uim3 + (7./60.)*uim2 + (-23./60.)*uim1 + (19./20.)*ui + (11./30.)*uip1 + (-1./30.)*uip2;  
        u4 = (1./60.)*uim2 + (-2./15.)*uim1 + (37./60.)*ui   + (37./60.)*uip1 + (-2./15.)*uip2 + (1./60.)*uip3;  
        u5 = (-1./30.)*uim1 + (11./30.)*ui  + (19./20.)*uip1 + (-23./60.)*uip2 + (7./60.)*uip3 + (-1./60.)*uip4;  
        u6 = (1./6.)*ui   + (29./20.)*uip1 + (-21./20.)*uip2 + (37./60.)*uip3 + (-13./60.)*uip4 + (1./30.)*uip5;  

        omt0 = g0*pow(eps+beta0, -power_param);
        omt1 = g1*pow(eps+beta1, -power_param);
        omt2 = g2*pow(eps+beta2, -power_param);
        omt3 = g3*pow(eps+beta3, -power_param);
        omt4 = g4*pow(eps+beta4, -power_param);
        omt5 = g5*pow(eps+beta5, -power_param);
        omts = omt0 + omt1 + omt2 + omt3 + omt4 + omt5;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4 + omt4*u5 + omt5*u6 )/omts);
    }

}

// ---------------------------------------------------- // 
// SECTION: WENO-Z reconstruction                       //
// ---------------------------------------------------- // 

// See: D.S. Balsara and C.-W. Shu, "Monotonicity Preserving Weighted 
//      Essentially Non-oscillatory Schemes with Increasingly High Order of 
//      Accuracy". J. Comput. Phys. 160 (2000), pp. 405-452.
// and other references for examples.
void WenoReconstruct_Z5( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq( g.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    double g0, g1, g2;           
    g0 = 0.1; g1 = 0.6; g2 = 0.3;

    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    const int meqn = g.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // Compute smoothness indicators (identical for left/right values):
        beta0 =(13./12.)*pow(uim2-2*uim1+ui,2)+0.25*pow(uim2-4*uim1+3*ui,2);
        beta1 =(13./12.)*pow(uim1-2*ui+uip1,2)+0.25*pow(uim1-uip1,2);
        beta2 =(13./12.)*pow(ui-2*uip1+uip2,2)+0.25*pow(3*ui-4*uip1+uip2,2);

        double tau5 = fabs( beta0 - beta2 );

        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
        u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
        u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;

        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*(1. + pow( tau5 / (eps+beta0), power_param ) );
        omt1 = g1*(1. + pow( tau5 / (eps+beta1), power_param ) );
        omt2 = g2*(1. + pow( tau5 / (eps+beta2), power_param ) );
        omts = omt0+omt1+omt2;

        // # Return 5th-order conservative reconstruction
        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3)/omts );

    }


}

void WenoReconstruct_Z7( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 7);
    double uim3, uim2, uim1, ui, uip1, uip2, uip3;
    double u1, u2, u3, u4;
    
    double beta0, beta1, beta2, beta3;
    double omt0, omt1, omt2, omt3, omts;

    const double g0 = 1./35., g1 = 12./35., g2 = 18./35., g3 = 4./35.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++){
        uim3 = g.get(m, 1);
        uim2 = g.get(m, 2);
        uim1 = g.get(m, 3);
        ui   = g.get(m, 4);
        uip1 = g.get(m, 5);
        uip2 = g.get(m, 6);
        uip3 = g.get(m, 7);

        beta0 = uim3*(  547.*uim3 -  3882.*uim2  + 4642.*uim1 - 1854.*ui) + 
                uim2*( 7043.*uim2 - 17246.*uim1  + 7042.*ui) +              
                uim1*(11003.*uim1 -  9402.*ui  ) + 2107. * pow(ui, 2);
        beta1 = uim2*(  267.*uim2 - 1642.*uim1   + 1602.*ui - 494.*uip1) + 
                uim1*( 2843.*uim1 - 5966.*ui     + 1922.*uip1) +
                ui*(   3443.*ui   - 2522.*uip1 ) + 547. * pow(uip1, 2);
        beta2 = uim1*( 547.*uim1 - 2522.*ui     + 1922.*uip1 - 494.*uip2) + 
                ui  *(3443.*ui   - 5966.*uip1   + 1602.*uip2 )   + 
                uip1*(2843.*uip1 - 1642.*uip2 ) + 267.*pow(uip2, 2);
        beta3 = ui*  ( 2107.*ui   -  9402.*uip1   + 7042.*uip2 - 1854.*uip3 ) + 
                uip1*(11003.*uip1 - 17246.*uip2   + 4642.*uip3 )              + 
                uip2*( 7043.*uip2 -  3882.*uip3 ) + 547.*pow(uip3, 2);

        double tau5 = abs( beta0 + 3.*(beta1-beta2) - beta3 );

        u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui;
        u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1;
        u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2;
        u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3;

        omt0 = g0*(1. + pow( tau5 / (eps+beta0), power_param ) );
        omt1 = g1*(1. + pow( tau5 / (eps+beta1), power_param ) );
        omt2 = g2*(1. + pow( tau5 / (eps+beta2), power_param ) );
        omt3 = g3*(1. + pow( tau5 / (eps+beta3), power_param ) );
        omts = omt0 + omt1 + omt2 + omt3;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4)/omts);
    }

}

void WenoReconstruct_Z9( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq(g.getsize(2), 9);
    double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4;
    double u1, u2, u3, u4, u5;
    
    double beta0, beta1, beta2, beta3, beta4;
    double omt0, omt1, omt2, omt3, omt4, omts;

    const double g0 = 1.0/126.0, g1 = 10./63., g2 = 10./21., g3 = 20./63., g4 = 5./126.;

    const int meqn = g.getsize(1);
    
    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    for(int m = 1; m <= meqn; m++)
    {

        uim4 = g.get(m, 1);
        uim3 = g.get(m, 2);
        uim2 = g.get(m, 3);
        uim1 = g.get(m, 4);
        ui   = g.get(m, 5);
        uip1 = g.get(m, 6);
        uip2 = g.get(m, 7);
        uip3 = g.get(m, 8);
        uip4 = g.get(m, 9);

        beta0 = uim4*(22658.*uim4 - 208501.*uim3 + 364863.*uim2 - 288007.*uim1 + 86329.*ui) +
                uim3*(482963.*uim3 - 1704396.*uim2 + 1358458.*uim1 - 411487. * ui) +
                uim2*(1521393.*uim2 - 2462076.*uim1 + 758823.*ui) + 
                uim1*(1020563.*uim1 - 649501.*ui) +
                107918.*pow(ui, 2);
        beta1 = uim3*(6908.*uim3  - 60871.*uim2 + 99213.*uim1 - 70237.*ui + 18079.*uip1) +
                uim2*(138563.*uim2 - 464976.*uim1 + 337018.*ui - 88297.*uip1) +
                uim1*(406293.*uim1 - 611976*ui + 165153.*uip1) +
                ui*(242723.*ui - 140251.*uip1) + 
                22658.*pow(uip1, 2);
        beta2 = uim2*(6908.*uim2 - 51001.*uim1 + 67923.*ui - 38947.*uip1 + 8209.*uip2) +
                uim1*(104963.*uim1 - 299076.*ui + 179098.*uip1 - 38947.*uip2) +
                ui*(231153.*ui - 299076.*uip1 + 67923.*uip2) +
                uip1*(104963.*uip1 - 51001.*uip2) +
                6908.*pow(uip2, 2);
        beta3 = uim1*(22658.*uim1 - 140251.*ui + 165153.*uip1 - 88297.*uip2 + 18079.*uip3) +
                ui*(242723.*ui - 611976.*uip1 + 337018.*uip2 - 70237.*uip3) +
                uip1*(406293.*uip1 - 464976.*uip2 + 99213.*uip3) +
                uip2*(138563.*uip2 - 60871.*uip3) +
                6908.*pow(uip3, 2);
        beta4 = ui*(107918.*ui - 649501.*uip1 + 758823.*uip2 - 411487.*uip3 + 86329.*uip4) +
                uip1*(1020563.*uip1 - 2462076.*uip2 + 1358458.*uip3 - 288007*uip4) + 
                uip2*(1521393.*uip2 - 1704396.*uip3 + 364863.*uip4) +
                uip3*(482963.*uip3 - 208501.*uip4) +
                22658.*pow(uip4, 2);

        // TODO - check this!
        double tau = abs( beta4 - beta0 );

        u1 = (1./5.)*uim4   + (-21./20.)*uim3 + (137./60.)*uim2 + (-163./60.)*uim1 + (137./60.)*ui;
        u2 = (-1./20.)*uim3 + (17./60.)*uim2  + (-43./60.)*uim1 + (77./60.)*ui     + (1./5.)*uip1;
        u3 = (1./30.)*uim2  + (-13./60.)*uim1 + (47./60.)*ui    + (9./20.)*uip1    + (-1./20.)*uip2;
        u4 = (-1./20.)*uim1 + (9./20.)*ui     + (47./60.)*uip1  + (-13./60.)*uip2  + (1./30.)*uip3;
        u5 = (1./5.)*ui     + (77./60.)*uip1  + (-43./60.)*uip2 + (17./60.)*uip3   + (-1./20.)*uip4;

        omt0 = g0*(1. + pow( tau / (eps+beta0), power_param ) );
        omt1 = g1*(1. + pow( tau / (eps+beta1), power_param ) );
        omt2 = g2*(1. + pow( tau / (eps+beta2), power_param ) );
        omt3 = g3*(1. + pow( tau / (eps+beta3), power_param ) );
        omt4 = g4*(1. + pow( tau / (eps+beta4), power_param ) );
        omts = omt0 + omt1 + omt2 + omt3 + omt4;

        g_reconst.set(m, 1, (omt0*u1 + omt1*u2 + omt2*u3 + omt3*u4 + omt4*u5)/omts);
    }

}

// ---------------------------------------------------- // 
// SECTION: Upwinded reconstruction with linear weights //
// ---------------------------------------------------- // 

// Conservative reconstruction based on the *linear* weights.  Do not use for
// problems with shocks.
void WenoReconstruct_FD5( const dTensor2& g, dTensor2& g_reconst )
{

    assert_eq( g.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    // double g0, g1, g2;           
    // g0 = 0.1; g1 = 0.6; g2 = 0.3;

    const int meqn = g.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // -- central finite difference reconstruction -- //
        u1 = ( 1./3.)*uim2 - (7./6.)*uim1 + (11./6.)*ui;
        u2 = (-1./6.)*uim1 + (5./6.)*ui   + ( 1./3.)*uip1;
        u3 = ( 1./3.)*ui   + (5./6.)*uip1 - ( 1./6.)*uip2;

        // 5th-order reconstruction (linear weights)
        g_reconst.set(m, 1, 0.1*u1+0.6*u2+0.3*u3 );

    }

}

void WenoReconstruct_FD7( const dTensor2& g, dTensor2& g_reconst )
{
    assert_eq(g.getsize(2), 7);
    double uim3, uim2, uim1, ui, uip1, uip2, uip3;
    double u1, u2, u3, u4;
    
    double beta0, beta1, beta2, beta3;
    double omt0, omt1, omt2, omt3, omts;

    const double g0 = 1./35., g1 = 12./35., g2 = 18./35., g3 = 4./35.;

    const int meqn = g.getsize(1);
    
    for(int m = 1; m <= meqn; m++)
    {

        uim3 = g.get(m, 1);
        uim2 = g.get(m, 2);
        uim1 = g.get(m, 3);
        ui   = g.get(m, 4);
        uip1 = g.get(m, 5);
        uip2 = g.get(m, 6);
        uip3 = g.get(m, 7);

        u1 = (-1./4. )*uim3 + (13./12.)*uim2 - (23./12.)*uim1 + (25./12.)*ui;
        u2 = ( 1./12.)*uim2 - ( 5./12.)*uim1 + (13./12.)*ui   + ( 1./4. )*uip1;
        u3 = (-1./12.)*uim1 + ( 7./12.)*ui   + ( 7./12.)*uip1 - ( 1./12.)*uip2;
        u4 = ( 1./4. )*ui   + (13./12.)*uip1 - ( 5./12.)*uip2 + ( 1./12.)*uip3;

        g_reconst.set(m, 1, (g0*u1 + g1*u2 + g2*u3 + g3*u4) );
   }
}

// 9th-order reconstruction (linear weights)
void WenoReconstruct_FD9( const dTensor2& g, dTensor2& g_reconst )
{
    assert_eq(g.getsize(2), 9);
    double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4;
    double u1, u2, u3, u4, u5;
    
    double beta0, beta1, beta2, beta3, beta4;
    double omt0, omt1, omt2, omt3, omt4, omts;

    const double g0 = 1.0/126.0, g1 = 10./63., g2 = 10./21., g3 = 20./63., g4 = 5./126.;

    const int meqn = g.getsize(1);
    for(int m = 1; m <= meqn; m++)
    {
        uim4 = g.get(m, 1);
        uim3 = g.get(m, 2);
        uim2 = g.get(m, 3);
        uim1 = g.get(m, 4);
        ui   = g.get(m, 5);
        uip1 = g.get(m, 6);
        uip2 = g.get(m, 7);
        uip3 = g.get(m, 8);
        uip4 = g.get(m, 9);

        u1 = (1./5.)*uim4   + (-21./20.)*uim3 + (137./60.)*uim2 + (-163./60.)*uim1 + (137./60.)*ui;
        u2 = (-1./20.)*uim3 + (17./60.)*uim2  + (-43./60.)*uim1 + (77./60.)*ui     + (1./5.)*uip1;
        u3 = (1./30.)*uim2  + (-13./60.)*uim1 + (47./60.)*ui    + (9./20.)*uip1    + (-1./20.)*uip2;
        u4 = (-1./20.)*uim1 + (9./20.)*ui     + (47./60.)*uip1  + (-13./60.)*uip2  + (1./30.)*uip3;
        u5 = (1./5.)*ui     + (77./60.)*uip1  + (-43./60.)*uip2 + (17./60.)*uip3   + (-1./20.)*uip4;

        g_reconst.set(m, 1, (g0*u1 + g1*u2 + g2*u3 + g3*u4 + g4*u5) );
    }
}

void WenoReconstructLLF( const dTensor2& g, dTensor2& g_reconst )
{
    
    const int meqn = g.getsize(1);
    for(int m = 1; m <= meqn; m++)
    {
        g_reconst.set(m, 1, g.get(m,1) );
    }

}

// ------------------------------------------------- // 
// SECTION: Central Finite difference approximations //
// ------------------------------------------------- // 

// First-derivative ( using a 5 point central stencil -> fourth-order )
void Diff1( double dx, const dTensor2& f, dTensor1& fx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils:
    assert_eq( f.getsize( 2 ), 5 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = (  f.get( m, 1 ) - f.get( m, 5 ) )*(1.0/12.0);
        tmp       += ( -f.get( m, 2 ) + f.get( m, 4 ) )*(2.0/ 3.0);
        fx.set( m, tmp / dx );
    }

}

// Central Finite difference approximations:
//
// First-derivative (using a 5 point central stencil)
//
// This is the scalar version of the above routine
double Diff1( double dx, 
    double f1, double f2, double f3, double f4, double f5 )
{

    double tmp = (  f1 - f5 )*(1.0/12.0);
    tmp       += ( -f2 + f4 )*(2.0/ 3.0);
    return tmp/dx;

}

// First-derivative non-conservative, based on WENO differentation (not
// WENO-reconstruction).
void Diff1NC( double dx, const dTensor2& g, dTensor1& fx )
{

    assert_eq( g.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    double g0, g1, g2;           
    g0 = 1./6.; g1 = 2./3.; g2 = 1./6.;

    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    const int meqn = g.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = g.get(m,1);
        uim1 = g.get(m,2);
        ui   = g.get(m,3);
        uip1 = g.get(m,4);
        uip2 = g.get(m,5);

        // Compute smoothness indicators (identical for left/right values):
        beta0 = pow(uim2-2.*uim1+ui,2);
        beta1 = pow(uim1-2.*ui+uip1,2);
        beta2 = pow(ui-2.*uip1+uip2,2);

        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 0.5  )*uim2 - (2.   )*uim1 + ( 1.5  )*ui;
        u2 = (-0.5  )*uim1 + (0.   )*ui   + ( 0.5  )*uip1;
        u3 = (-1.5  )*ui   + (2.   )*uip1 - ( 0.5  )*uip2;
        
        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*pow(eps+beta0,-power_param);
        omt1 = g1*pow(eps+beta1,-power_param);
        omt2 = g2*pow(eps+beta2,-power_param);
        omts = omt0+omt1+omt2;

        // # Return 5th-order approximation to derivative
        fx.set(m, ( (omt0*u1 + omt1*u2 + omt2*u3)/omts) / dx );

    }

}


// Central Finite difference approximations:
//
// Second-derivative (using a 5 point central stencil)
void Diff2( double dx, const dTensor2& f, dTensor1& fxx )
{

    const int mcomps = f.getsize( 1 );

    // TODO - include options for larger stencils
    assert_eq( f.getsize( 2 ), 5 );
    for( int m=1; m <= mcomps; m++ )
    {
        double tmp = ( -f.get( m, 1 ) - f.get( m, 5 ) )*(1.0/12.0);
        tmp       += (  f.get( m, 2 ) + f.get( m, 4 ) )*(4.0/ 3.0);
        tmp       += ( -f.get( m, 3 )                 )*(5.0/ 2.0);
        fxx.set( m, tmp / (dx*dx) );
    }

}

// First-derivative non-conservative, based on WENO differentation (not
// WENO-reconstruction).
void Diff2NC( double dx, const dTensor2& f, dTensor1& fxx )
{

    assert_eq( f.getsize(2), 5 );

    // Stencil and the smaller three point derivatives:
    double uim2,uim1,ui,uip1,uip2;
    double u1,u2,u3;

    double beta0, beta1, beta2;  // smoothness indicators
    double omt0, omt1, omt2, omts;

    // linear weights
    const double g0 = -1./12.; 
    const double g1 =  7./6.; 
    const double g2 = -1./12.;

    const double eps         = global_ini_params.get_epsilon();       
    const double power_param = global_ini_params.get_power_param();   // Default: p=2

    const int meqn = f.getsize(1);
    for( int m=1; m <= meqn; m++ )
    {

        uim2 = f.get(m,1);
        uim1 = f.get(m,2);
        ui   = f.get(m,3);
        uip1 = f.get(m,4);
        uip2 = f.get(m,5);

        // Compute smoothness indicators (identical for left/right values):
        beta0 = pow(uim2-2.*uim1+ui,2);
        beta1 = pow(uim1-2.*ui+uip1,2);
        beta2 = pow(ui-2.*uip1+uip2,2);

        // 3rd-order reconstructions using small 3-point stencils
        u1 = ( 1.0  )*uim2 - (2.   )*uim1 + ( 1.0  )*ui;
        u2 = ( 1.0  )*uim1 - (2.   )*ui   + ( 1.0  )*uip1;
        u3 = ( 1.0  )*ui   - (2.   )*uip1 + ( 1.0  )*uip2;
        
        // Compute nonlinear weights and normalize their sum to 1
        omt0 = g0*pow(eps+beta0,-power_param);
        omt1 = g1*pow(eps+beta1,-power_param);
        omt2 = g2*pow(eps+beta2,-power_param);
        omts = omt0+omt1+omt2;

        // # Return 5th-order approximation to derivative
        fxx.set(m, ( (omt0*u1 + omt1*u2 + omt2*u3)/omts) / (dx*dx) );

    }

}
#ifndef _WENORECONSTRUCT_H_
#define _WENORECONSTRUCT_H_

// -- Jiang and Shu reconstructions -- //
void WenoReconstruct_JS5 ( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_JS7 ( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_JS9 ( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_JS11( const dTensor2& g, dTensor2& g_reconst );

// -- WENO-Z reconstructions -- //
void WenoReconstruct_Z5( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_Z7( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_Z9( const dTensor2& g, dTensor2& g_reconst );

// -- Reconstructions based on linear weights -- //
// These are useful for running convergence studies, and should not be used
// for problems with discontinuities.
void WenoReconstruct_FD5( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_FD7( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_FD9( const dTensor2& g, dTensor2& g_reconst );

// -- Lax-Friedrichs solver -- //
// Although this is not a "WENO" reconstruction, we include this here so that
// we can use the same call to ConstructL for all solvers
void WenoReconstructLLF( const dTensor2& g, dTensor2& g_reconst );

// Wrapper function that provides access to each of the above through looking
// at the global variable wenoParams.
//void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
typedef void (*reconstruct_t)(const dTensor2&, dTensor2&);
reconstruct_t GetWenoReconstruct();

#endif
#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

const double pi    = 3.1415926535897932385; // M_PI; // 4.0*atan(1.0);

// Common square roots encountered in the code:
const double sq2   = 1.4142135623730950488;
const double sq3   = 1.7320508075688772935;
const double sq5   = 2.2360679774997896964;
const double sq6   = 2.4494897427831780982;
const double sq7   = 2.6457513110645905905;
const double sq10  = 3.1622776601683793320;
const double sq11  = 3.3166247903553998491;
const double sq13  = 3.6055512754639892931;
const double sq15  = 3.8729833462074168852;
const double sq19  = 4.3588989435406735522;
const double sq23  = 4.7958315233127195416;
const double sq30  = 5.4772255750516611346;
const double sq35  = 5.9160797830996160426;
const double sq71  = 8.4261497731763586306;
const double sq105 = 1.0246950765959598383e+01;

// Common recipricols of the above square roots enountered in the code:
const double osq2  = 7.0710678118654752440e-01;
const double osq3  = 5.7735026918962576449e-01;
const double osq5  = 4.4721359549995793928e-01;
const double osq7  = 3.7796447300922722721e-01;
const double osq13 = 2.7735009811261456101e-01;
const double osq19 = 2.2941573387056176591e-01;

// Common fractions where we would prefer multiplication over division:
const double onehalf    = 0.5;
const double onethird   = 3.3333333333333333333e-01;
const double oneseventh = 1.4285714285714285714e-01;
const double oneninth   = 1.1111111111111111111e-01;

// Used to test for machine zero:
const double EPSILON = 0.5e-15;

// Faster lookup for small factorials.  See also: dog_math.h
const double  factorial[11] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0 };
const int    ifactorial[11] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 };

#endif
#include "dog_io.h"
#include "debug.h"
#define NDIMS 0
#include "tensors1d.h"
#include <unistd.h> // for unlink
#include <string>

// Is there a higher-level system command to do a file copy?
void copyFile(const char* inFileName, const char* outFileName)
{
  FILE *inFile = fopen(inFileName,"r");
  if(!inFile) eprintf("could not open file for reading: %s",inFileName);
  FILE *outFile = fopen(outFileName,"w");
  if(!outFile) eprintf("could not open file for writing: %s",outFileName);
  int c;
  while(c = fgetc(inFile), c!=EOF)
  {
    fputc(c,outFile);
  }
  fclose(inFile);
  fclose(outFile);
}

// returns true if file exists
//
// need to change this to call opendir so that the directory
// cache will be updated. otherwise this doesn't work over a
// networked file system, if matlab (the requesting process)
// is running on a different machine.  But then we should
// really be using tcpip to transfer data anyway rather than
// writing and reading files.
//
bool flagfile_check(const char* filename)
{
  //DIR *dp = opendir (get_outputdir());
  FILE* file = fopen(filename,"r");
  if(!file) return false;
  fclose(file);
  return true;
}

// returns 0 if successfully created flagfile
int flagfile_create(const char* filename)
{
  FILE* file = fopen(filename,"w");
  if(!file) return 1;
  fclose(file);
  return 0;
}

// return 0 if successfully removed flagfile
int flagfile_remove(const char* filename)
{
  int err = unlink(filename);
  return err;
}

#if 0
void create_lock_file(const char* basefilename)
{
  using namespace std;
  string lockfilename = string(basefilename)+".lock";
  //int pfd = open(lockfilename.c_str(),
  //  O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  FILE* file = fopen(lockfilename.c_str(),"w");
  fprintf(file,"This lock file should not persist; you can delete it.\n");
  if(!file) eprintf("Could not create lock file %s", lockfilename.c_str());
  fclose(file);
}

void remove_lock_file(const char* basefilename)
{
  using namespace std;
  string lockfilename = string(basefilename)+".lock";
  int error = unlink(lockfilename.c_str());
  if(error) eprintf("could not unlink file %s", lockfilename.c_str());
}
#endif

void fprint_array(FILE* file, const int* arr, int first_idx, int last_idx,
  char field_sep, const char*end)
{
  if(last_idx-first_idx>=0) fprintf(file,"%d",arr[first_idx]);
  for(int i=first_idx+1;i<=last_idx;i++)
  {
    fputc(field_sep,file);
    fprintf(file,"%d",arr[i]);
  }
  fprintf(file,"%s",end);
}

void fprint_tensor(FILE* file, const iTensorBase& t,
  char field_sep, const char*end)
{
  if(t.numel()>0) fprintf(file,"%d",t.vget(0));
  for(int i=1;i<t.numel();i++)
  {
    fputc(field_sep,file);
    fprintf(file,"%d",t.vget(i));
  }
  fprintf(file,"%s",end);
}

void fprint_tensor(FILE* file, const dTensorBase& t,
  char field_sep, const char*end)
{
  if(t.numel()>0) fprintf(file,"%24.16e",t.vget(0));
  for(int i=0;i<t.numel();i++)
  {
    fputc(field_sep,file);
    fprintf(file,"%24.16e",t.vget(i));
  }
  fprintf(file,"%s",end);
}
#ifndef __DOG_IO_H__
#define __DOG_IO_H__

#include <stdio.h>

// Is there a higher-level system command to do a file copy?
void copyFile(const char* inFileName, const char* outFileName);

// conversion between strings and arrays
//
class iTensorBase;
class dTensorBase;
void fprint_tensor(FILE* file, const iTensorBase& t,
  char field_sep=',', const char*end="\n");
void fprint_tensor(FILE* file, const dTensorBase& t,
  char field_sep=',', const char*end="\n");
void fprint_array(FILE* file, const int* arr, int first_idx, int last_idx,
  char field_sep=',', const char*end="\n");

bool flagfile_check(const char* filename);
int flagfile_create(const char* filename);
int flagfile_remove(const char* filename);

#endif
#include <cmath>
#include "dog_math.h"
#define NDIMS 0
#include "dogdefs.h"

//////////////////////////////////////////////////////////////////////////////
//  The purpose of this module is to describe common math functions
//  what we'd like to be able to perform in DogPack.
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Returns the maximum value of a single vector
//////////////////////////////////////////////////////////////////////////////
double Max(const dTensor1& vec)
{
    double max = 0.0;
    for( int i=1; i <= vec.getsize(); i++)
    {
        max = Max(max, vec.get(i));
    }
    return max;
}

//////////////////////////////////////////////////////////////////////////////
// Returns the maximum value of a single vector
//////////////////////////////////////////////////////////////////////////////
/*double Max(const dTensorBC1& vec)
{
    double Max(double,double);
    double max = 0.0;
    double mbc = vec.getmbc();
    for( int i=1-mbc; i <= vec.getsize()+mbc; i++)
    {
        max = Max(max, vec.get(i));
    }
    return max;
}
*/

//////////////////////////////////////////////////////////////////////////////
// Returns a vector containing the maximum value of each component in vec1 and
// vec2
//////////////////////////////////////////////////////////////////////////////
void Max(const dTensor1& vec1, const dTensor1& vec2, dTensor1& max)
{
    assert_eq(vec1.getsize(),vec2.getsize());
    assert_eq(max.getsize(),vec1.getsize());

    for( int i=1; i <= vec1.getsize(); i++)
    {
        max.set(i,  Max(vec1.get(i), vec2.get(i)) );
    }
}

//////////////////////////////////////////////////////////////////////////////
// Returns the maximum value found in Matrix A
//////////////////////////////////////////////////////////////////////////////
double Max(const dTensorBC2& A)
{
    int mbc = A.getmbc();
    double tmp = A.get(1,1);
    for( int i=1-mbc; i <= A.getsize(1)+mbc; i++)
    for( int j=1-mbc; j <= A.getsize(2)+mbc; j++)
    {
        tmp = Max(tmp, A.get(i,j) );
    }
    return tmp;
}

// Compute the product: A*B. Save output in C.
//
// Tensors need to have the same size.  That is, if 
// A is mxn and B is nxl, then C should be mxl.
//
void TensorMultiply( const dTensor2& A, const dTensor2& B, dTensor2& C )
{
    // Quick error check:
    assert_eq( A.getsize( 2 ), B.getsize(1) );
    assert_eq( A.getsize( 1 ), C.getsize(1) );
    assert_eq( B.getsize( 2 ), C.getsize(2) );

    for( int k=1; k <= A.getsize(1); k++ )
    for( int l=1; l <= B.getsize(2); l++ )
    {

        // TODO: FOR LARGE MATRICES, IT MAY BE FASTER TO FIRST TAKE THE
        // TRANSPOSE OF B, AND THEN DO THE MULTIPLCIATION ... (-DS)
        double tmp = 0.;
        for( int j=1; j <= A.getsize(2); j++ )
        {
            tmp += A.get(k,j)*B.get(j, l );
        }
        
    }
}

// Compute the product: A*v. Save output in b.
//
void MatVecMultiply( const dTensor2& A, const dTensor1& v, dTensor1& b )
{
    // Quick error check:
    assert_eq( A.getsize( 2 ), v.getsize() );
    assert_eq( A.getsize( 1 ), b.getsize() );

    for( int k=1; k <= A.getsize(1); k++ )
    {

        double tmp = 0.;
        for( int j=1; j <= A.getsize(2); j++ )
        {
            tmp += A.get(k,j)*v.get(j);
        }
        b.set( k, tmp );
        
    }
}

//////////////////////////////////////////////////////////////////////////////
// Multiply two functions Q1 * Q2 and write to output qnew
//
// Q1(x) = \sum_{i,k} Q1_i^{(k)} \phi_i^{(k)}(x)
// Q2(x) = \sum_{i,k} Q2_i^{(k)} \phi_i^{(k)}(x)
//
// The resultant function qnew is after doing the projections.
//////////////////////////////////////////////////////////////////////////////
void MultiplyFunctions(const dTensorBC3& Q1, 
    const dTensorBC3& Q2, dTensorBC3& qout)
{
    const int melems = qout.getsize(1);
    const int meqn   = qout.getsize(2);
    const int kmax   = qout.getsize(3);
    const int mbc    = qout.getmbc();

    // TODO - I THINK THIS COULD BE MUCH MORE EFFICIENT IF WE WROTE THESE OUT BY
    // HAND RATHER THAN USING A MATRIX TO DO THE MULTIPLICATION.  THERE ARE A
    // LOT OF ZEROS IN THIS MATRIX THAT DO NOT NEED TO BE INCLUDED
    dTensor3 M(5,5,5);
    M.setall(0.);

    // k == 5
    M.set(5,1, 5, 1.0 );
    M.set(1,5, 5, 1.0 );

    M.set(2,4, 5, 4.0*sq3*sq7/21.0 );
    M.set(4,2, 5, M.get(2,4,5) );

    M.set(3,5, 5, 20.0/77.0*sq5 );
    M.set(5,3, 5, M.get(3,5,5) );

    M.set(3,3, 5, 6.0/7.0 );
    M.set(4,4, 5, 6.0/11.0 );
    M.set(5,5, 5, 486.0/1001.0 );

    // k == 4
    M.set(1,4, 4, 1.0 );
    M.set(4,1, 4, 1.0 );
    M.set(2,3, 4, sq3*sq5*sq7*3.0/35.0 );
    M.set(3,2, 4, M.get(2,3,4) );
    M.set(2,5, 4, 4.0*sq3*sq7/21.0 );
    M.set(5,2, 4, 4.0*sq3*sq7/21.0 );
    M.set(3,4, 4, 4.0*sq5/15.0 );
    M.set(4,3, 4, 4.0*sq5/15.0 );
    M.set(4,5, 4, 6.0/11.0 );
    M.set(5,4, 4, 6.0/11.0 );

    // k == 3
    M.set(1,3, 3, 1.0 );
    M.set(3,1, 3, 1.0 );
    M.set(2,2, 3, 2.0*sq5/5.0 );
    M.set(2,4, 3, 3.0*sq3*sq5*sq7/35.0 );
    M.set(4,2, 3, 3.0*sq3*sq5*sq7/35.0 );
    M.set(3,3, 3, 2.0*sq5/7.0 );
    M.set(3,5, 3, 6.0/7.0 );
    M.set(5,3, 3, 6.0/7.0 );
    M.set(4,4, 3, 4.0*sq5/15.0 );
    M.set(5,5, 3, 20.0*sq5/77.0 );

    // k == 2
    M.set(1,2, 2, 1.0 );
    M.set(2,1, 2, 1.0 );
    M.set(2,3, 2, 2.0*sq5/5.0 );
    M.set(3,2, 2, 2.0*sq5/5.0 );
    M.set(3,4, 2, sq3*sq5*sq7*3.0/35.0 );
    M.set(4,3, 2, sq3*sq5*sq7*3.0/35.0 );
    M.set(4,5, 2, sq3*sq7*4.0/21.0 );
    M.set(5,4, 2, sq3*sq7*4.0/21.0 );

    // k == 1
    M.set(1,1, 1, 1.0 );
    M.set(2,2, 1, 1.0 );
    M.set(3,3, 1, 1.0 );
    M.set(4,4, 1, 1.0 );
    M.set(5,5, 1, 1.0 );

//#pragma omp parallel for
    for(int i=1-mbc; i<=melems+mbc; i++ )
    for(int me=1; me <= meqn; me++)
    {
        double tmp = 0.0;
        for( int  k=1; k <= kmax; k++)
        {
            double tmp = 0.0;
            for( int k1=1; k1 <= kmax; k1++)
            for( int k2=1; k2 <= kmax; k2++)
            { tmp += Q1.get(i, me, k1)*Q2.get(i, me, k2)*M.get(k1,k2,k); }
            qout.set(i,me,k,tmp);
        }

    }

}

/*
        switch(kmax)
        {

            case 5:

                ////////////// k == 5 /////////////
                tmp = 0.0;

                tmp += 486.0/1001.0*Q1.get(i,me,5)*Q2.get(i,me,5);
                tmp += 6.0/11.0*Q1.get(i,me,4)*Q2.get(i,me,4);
                tmp += 6.0/7.0*Q1.get(i,me,3)*Q2.get(i,me,3);

                tmp += Q1.get(i,me,5) * Q2.get(i,me,3) * 
                    20.0 * sq5 / 77.0;
                tmp += Q2.get(i,me,5) * Q1.get(i,me,3) * 
                    20.0 * sq5 / 77.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,4) * 
                    4.0*sq3*sq7/21.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,4) * 
                    4.0*sq3*sq7/21.0;

                tmp += Q1.get(i,me,5) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,5) * Q1.get(i,me,1);
                qnew.set(i,me,5,tmp);


                ////////////// k == 4 /////////////
                tmp = 0.0;
                tmp += Q1.get(i,me,2) * Q2.get(i,me,5) * 
                    4.0 * sq3 * sq7 / 21.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,5) * 
                    4.0 * sq3 * sq7 / 21.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,5) * 
                    6.0 / 11.0;
                tmp += Q2.get(i,me,4) * Q1.get(i,me,5) * 
                    6.0 / 11.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,3) * 
                    4.0 * sq5 / 15.0;
                tmp += Q2.get(i,me,4) * Q1.get(i,me,3) * 
                    4.0 * sq5 / 15.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,3) * 
                    3.0*sq3*sq5*sq7 / 35.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,3) * 
                    3.0*sq3*sq5*sq7 / 35.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,4) * Q1.get(i,me,1);
                qnew.set(i,me,4,tmp);


                ////////////// k == 3 /////////////
                tmp = 0.0;

                tmp += 20.0*sq5/77.0*Q1.get(i,me,5)*Q2.get(i,me,5);
                tmp += Q1.get(i,me,3) * Q2.get(i,me,5) * 
                    6.0 / 7.0;
                tmp += Q2.get(i,me,3) * Q1.get(i,me,5) * 
                    6.0 / 7.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,4) * 
                        4.0 / (3.0*sq5);

                tmp += 2.0*sq5/7.0*Q1.get(i,me,3)*Q2.get(i,me,3);

                tmp += Q1.get(i,me,2) * Q2.get(i,me,4) * 
                    3.0*sq5*sq3*sq7 / 35.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,4) * 
                    3.0*sq5*sq3*sq7 / 35.0;

                tmp += 2.0 / sq5 * Q1.get(i,me,2)*Q2.get(i,me,2);

                tmp += Q1.get(i,me,3) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,3) * Q1.get(i,me,1);
                qnew.set(i,me,3,tmp);

                ////////////// k == 2 /////////////
                tmp = 0.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,5) * 4.0*sq3*sq7/21.0;
                tmp += Q2.get(i,me,4) * Q1.get(i,me,5) * 4.0*sq3*sq7/21.0;

                tmp += Q1.get(i,me,3) * Q2.get(i,me,4) * 
                    3.0*sq3*sq5*sq7/35.0;
                tmp += Q2.get(i,me,3) * Q1.get(i,me,4) * 
                    3.0*sq3*sq5*sq7/35.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,3) * 2.0 / sq5;
                tmp += Q1.get(i,me,3) * Q2.get(i,me,2) * 2.0 / sq5;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,1);
                tmp += Q1.get(i,me,1) * Q2.get(i,me,2);
                qnew.set(i,me,2,tmp);

                ////////////// k == 1 /////////////
                tmp = 0.0;
                for( int k=1; k <= kmax; k++ )
                { tmp += Q1.get(i,me,k) * Q2.get(i,me,k); }
                qnew.set(i,me,1, tmp);
                break;


            case 4:

                ////////////// k == 4 /////////////
                tmp = 0.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,3) * 
                    4.0 * sq5 / 15.0;
                tmp += Q2.get(i,me,4) * Q1.get(i,me,3) * 
                    4.0 * sq5 / 15.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,3) * 
                    3.0*sq3*sq5*sq7 / 35.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,3) * 
                    3.0*sq3*sq5*sq7 / 35.0;

                tmp += Q1.get(i,me,4) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,4) * Q1.get(i,me,1);
                qnew.set(i,me,4,tmp);


                ////////////// k == 3 /////////////
                tmp = 0.0;
                tmp += Q1.get(i,me,4) * Q2.get(i,me,4) * 
                        4.0 / (3.0*sq5);

                tmp += 2.0*sq5/7.0*Q1.get(i,me,3)*Q2.get(i,me,3);

                tmp += Q1.get(i,me,2) * Q2.get(i,me,4) * 
                    3.0*sq5*sq3*sq7 / 35.0;
                tmp += Q2.get(i,me,2) * Q1.get(i,me,4) * 
                    3.0*sq5*sq3*sq7 / 35.0;

                tmp += 2.0 / sq5 * Q1.get(i,me,2)*Q2.get(i,me,2);

                tmp += Q1.get(i,me,3) * Q2.get(i,me,1);
                tmp += Q2.get(i,me,3) * Q1.get(i,me,1);
                qnew.set(i,me,3,tmp);

                ////////////// k == 2 /////////////
                tmp = 0.0;

                tmp += Q1.get(i,me,3) * Q2.get(i,me,4) * 
                    3.0*sq3*sq5*sq7/35.0;
                tmp += Q2.get(i,me,3) * Q1.get(i,me,4) * 
                    3.0*sq3*sq5*sq7/35.0;

                tmp += Q1.get(i,me,2) * Q2.get(i,me,3) * 2.0 / sq5;
                tmp += Q1.get(i,me,3) * Q2.get(i,me,2) * 2.0 / sq5;
                tmp += Q1.get(i,me,2) * Q2.get(i,me,1);
                tmp += Q1.get(i,me,1) * Q2.get(i,me,2);
                qnew.set(i,me,2,tmp);

                ////////////// k == 1 /////////////
                tmp = 0.0;
                for( int k=1; k <= kmax; k++ )
                {
                    tmp += Q1.get(i,me,k) * Q2.get(i,me,k);
                }
                qnew.set(i,me,1, tmp);
                break;

            default:
                perror(" method not implemented yet\n");
                exit(1);
        }
    }
    */
#ifndef DOG_MATH_H
#define DOG_MATH_H

#include <cmath> // This header extends math.h
#include "constants.h"

// function declarations
// (only declaring the functions from dog_math.cpp that actually get used)

//////////////////////
// inline functions //
//////////////////////

// sgnum function.  undefined behavior for x near zero
inline double sgn( double x )
{ return x<0?-1.:1.; }

inline double Max(double a, double b)
{ return a>b?a:b; }

inline double Min(double a, double b)
{ return a<b?a:b; }

inline int iMax(int a, int b)
{ return a>b ? a : b; }

inline int iMin(int a, int b)
{ return a<b ? a : b; }

// Cleaner mod function.  
// Returns n % m, where the result is between 0 and m-1.
inline int iMod(int n, int m)
{
    int modval = n%m;
    return n%m < 0 ? (modval)+m : modval;
}

inline double minmod(double a, double b)
{
    if(std::signbit(a))
    {
      if(std::signbit(b)) return Max(a,b);
      return 0.;
    }
    if(std::signbit(b)) return 0.;
    return Min(a,b);
}

inline double minmod(double a, double b, double c)
{
  if(std::signbit(a))
  {
    if(std::signbit(b) && std::signbit(c)) return Max(Max(a,b),c);
    return 0;
  }
  if(std::signbit(b) || std::signbit(c)) return 0;
  return Min(Min(a,b),c);
}

// poor man's Factorial function //
// Note: if n < 0, this function returns 1.
//
// See: constants.h: factorial for a single access array of factorials
// 0! through 10!
inline int Factorial(int n)
{
    int result = 1;
    while( n > 1 )
    { result *= n--; }
    return result;
}

inline bool almost_eq( double a, double b )
{
    return fabs(a-b) < EPSILON ? true : false;
}
#endif // DOG_MATH_H
#include <stdio.h>
#include <cstring>     // for strcmp
#include "assert.h"
#include "dog_str.h"
#include "tensors1d.h"
#include "debug.h"

bool str_eq(const char* str1, const char* str2)
{
    // for case-insensitive matching change this to strcasecmp
    return !strcmp(str1,str2);
}

// convert comma-separated list to array of integers.
//
// user is responsible to
// delete [] int_array;
//
int new_array_from_str(
        const char* str_in,
        int* &int_array,
        int first_idx,
        char field_sep)
{
    assert(first_idx==0 || first_idx==1);

    int_array = 0;
    char *str = (char*)str_in;
    char *ptr = str;
    if(ptr==0) return -1;
    if(*ptr==0) return 0;

    // allocate int_array
    int number_of_numbers=0;
    {
        // count number of commas in string
        int number_of_commas=0;
        for(char *ptr=str; *ptr; ptr++) if(*ptr == field_sep) number_of_commas++;
        number_of_numbers = number_of_commas+1;
        int array_size = first_idx+number_of_numbers;
        int_array=new int[array_size];
        for(int i=0;i<array_size;i++) int_array[i]=0;
    }

    // populate int_array with numbers from list
    ptr=str;
    int idx=first_idx;
    while(*ptr) 
    {
        int val;
        int successful_match = sscanf(ptr,"%d",&val);
        if(!successful_match) {
            printf("problem: could not match %s as integer\n",ptr);
            delete [] int_array; int_array=0; return -1;
        } else {
            int_array[idx++]=val;
        }
        // advance ptr past next comma
        while(*ptr && *ptr++ != field_sep);
    }
    assert(idx==(first_idx+number_of_numbers));
    return number_of_numbers;
}

// int count_fields(const char* str_in, char field_sep)
// {
//   char *str = (char*)str_in;
//   char *ptr = str;
//   // a NULL or empty string is considered to have no fields
//   if(ptr==0) return 0;
//   if(*ptr==0) return 0;
//   int num_fields=1;
//   // count number of field_sep in string
//   for(char *ptr=str; *ptr; ptr++) if(*ptr == field_sep) num_fields++;
//   return num_fields;
// }
// 
// int count_char(const char* str, char c)
// {
//   if(str==0) return -1;
//   int num_c=0;
//   // count number of occurrences of c in string
//   for(char *ptr=(char*)str; *ptr; ptr++) if(*ptr == c) num_c++;
//   return num_c;
// }

// field_sep is typically ',' or '\n'
bool str_into_tensor(const char* str, iTensorBase& t, char field_sep)
{
    int idx=0;
    bool ret=true;
    char* ptr = (char*)str;
    while(*ptr && idx < t.numel()) 
    {
        int val;
        int successful_match = sscanf(ptr,"%d",&val);
        if(!successful_match) {
            Wprintf("could not match %s as integer\n",ptr);
            ret=false;
        } else {
            t.vset(idx++, val);
        }
        // advance ptr past next comma
        while(*ptr && *ptr++ != field_sep);
    }
    return ret;
};

// field_sep is typically ',' or '\n'
bool str_into_tensor(const char* str, dTensorBase& t, char field_sep)
{
    int idx=0;
    bool ret=true;
    char* ptr = (char*)str;
    while(*ptr && idx < t.numel()) 
    {
        double val;
        int successful_match = sscanf(ptr,"%lf",&val);
        if(!successful_match) {
            Wprintf("could not match %s as double\n",ptr);
            ret=false;
        } else {
            t.vset(idx++, val);
        }
        // advance ptr past next comma
        while(*ptr && *ptr++ != field_sep);
    }
    return ret;
};

#ifndef __DOG_STR_H__
#define __DOG_STR_H__

class iTensorBase;

// test strings for equality
bool str_eq(const char* str1, const char* str2);

// copy the contents of input string into the iTensor
bool str_into_tensor(const char* str, iTensorBase& t, char sep=',');

// Create a new array redirect the input pointer to this array from input string
// User is expected to delete the allocated memory when finished
int new_array_from_str( const char* str, int* &arr, int first_idx, char sep=',');

#endif
#ifndef _DOGDEFS_H_
#define _DOGDEFS_H_

// Handful of header files that are used in most files

#include "assert.h"        // includes <stdlib.h>
#include "debug.h"         // includes <stdio.h>
#include "constants.h"
#include "tensors.h"

#endif /*_DOGDEFS_H_*/
#include <string>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "dogdefs.h"
#include "dimdefs.h"
#include "IniParams.h"

/*
 * Common main function that's called by every application.
 *
 * Each application has its own main file, and Makefile that builds that local
 * main.  When the code is built in any application directory, every
 * application links to this common main function, that executes the code.
 *
 * The purpose of placing this extra layer between appname/main.cpp and
 * RunFinpack is to make the main function in each application as short as
 * possible.
 *
 */
int main_global(int argc, char* argv[])
{
    using std::string;
    using std::cout;
    using std::setprecision;
    using std::setw;
    using std::scientific;
    using std::endl;

    string parameters_ini_filename;

    parameters_ini_filename = 
        argc == 1 ? "parameters.ini" : argv[1];
       
    cout << "Running with configuration file: "
         << parameters_ini_filename
         << endl;

    global_ini_params.init(parameters_ini_filename);

    // Get current time
    double time1 = time(NULL);

    // Call startscript (Default: scripts/startscript, otherwise it checks for
    // a local file called 'startscript' from the application's directory)
    void RunStartScript(string parameters_ini_filename);
    RunStartScript(parameters_ini_filename);

    // Call the ``RunFinpack'' routine, which executes the code
    // Each dimension has its own version of this routine.
    int RunFinpack( );
    int m = RunFinpack( );

    // Get current time
    double time2 = time(NULL);

    // Output elapsed running time
    cout << setprecision(5);
    cout << " Total elapsed time in seconds = " << setw(11) 
        << scientific << time2-time1 << endl << endl;

    return m;

}
#include <string>
#include <iostream>
#include <stdlib.h>
#include "debug.h"
using namespace std;
//#include "stdio.h"

// so that we can print doubles to desired precision
//
void assert_error(const char* file, int line, const char* func,
  const char* op, const char* lhs_str, const char* rhs_str,
  double lhs, double rhs)
{
  fprintf(stderr,"ERROR in file %s, line %d, function %s"
      "\n\tassertion failed: %s %s %s, i.e., %24.16e %s %24.16e\n",
    file, line, func, lhs_str, op, rhs_str, lhs, op, rhs);
  abort();
}

#define implement_assert_errmsg(t1,t2) \
  void assert_error(const char* file, int line, const char* func, \
    const char* op, const char* lhs_str, const char* rhs_str, \
    t1 lhs, t2 rhs) \
  { \
    std::cerr<< "ERROR in file " << file << ", line " << line  \
      << ", function " << func  \
      <<"\n\tassertion failed: " << lhs_str << op << rhs_str \
      << ", i.e., " << lhs << op << rhs << endl; \
      abort(); \
  }

implement_assert_errmsg(int,int);
implement_assert_errmsg(const char*,const char*);
implement_assert_errmsg(const string&,const string&);
//implement_assert_errmsg(string,string);

//implement_assert_errmsg(double,double);
//void assert_error(const char* file, int line, const char* func, \
//  const char* op, const char* lhs_str, const char* rhs_str, \
//  double lhs, double rhs) \
//{ \
//  fprintf(stderr, "ERROR, %s(), %s:%d:\n\t"
//    "assertion failed: %s %s %s, i.e., %.16e %s %.16e\n",
//    func, file, line, 
//    lhs_str, op, rhs_str, lhs, op, rhs);
//    abort();
//}

// void assert_error(const char* file, int line, const char* func,
//   const char* op, const char* lhs_str, const char* rhs_str,
//   double lhs, double rhs)
// {
//   std::cerr<< "ERROR in file " << file << ", line " << line 
//     << ", function " << func 
//     <<"\n\tassertion failed: " << lhs_str << op << rhs_str
//     << ", i.e., " << lhs << op << rhs << endl;
//   abort();
// }

//  std::cerr<< "ERROR in file " << __FILE__ << ", line " << __LINE__ \
//           << ", function " << __func__ \
//           <<"\n\tassertion failed: " #lhs "==" #rhs ", lhs=" \
//           <<lhs<<", rhs="<<rhs<<std::endl; \

/*
 * This file is intended to supersede the system assert.h files
 * Unlike other header files, "assert.h" may usefully be included
 * multiple times, with and without NDEBUG defined.
 */

#include <stdio.h>          // for printf
#include <stdlib.h>         // for abort
//#include <sys/cdefs.h>

// asserts higher than MAX_ASSERT_LEVEL are compiled out
// (to avoid slowing execution)
// uncomment next line to make this variable global
// #undef MAX_ASSERT_LEVEL
#ifndef MAX_ASSERT_LEVEL
  #define MAX_ASSERT_LEVEL 2
#endif

// for consistency with functionality of assert.h
#ifdef NDEBUG // completely turns off assert statements
// if debug is turned off create vacuous versions
// of everything that user might use
// (a list of what we currently actually support;
// other stuff below is experimental)
#define assert(e) ((void)0)
#define assert1(e) ((void)0)
#define assert2(e) ((void)0)
#define assert3(e) ((void)0)
#define assert_printf(args...) ((void)0)
#define assert_printf1(args...) ((void)0)
#define assert_printf2(args...) ((void)0)
#define assert_printf3(args...) ((void)0)
#define assert_almost_eq(a) ((void)0)
#define assert_eq(a,b) ((void)0)
#define assert_ne(a,b) ((void)0)
#define assert_gt(a,b) ((void)0)
#define assert_lt(a,b) ((void)0)
#define assert_ge(a,b) ((void)0)
#define assert_le(a,b) ((void)0)
#define assert_isnum(a) ((void)0)

#else // ifndef NDEBUG

// override system assert.h
#define dassert_fileLine(e, file, line, func) \
    (void)(printf("ERROR: file %s, line %d, function %s:\n\tfailed assertion: (%s)\n", file, line, func,e),abort())
#define dassert_printf_fileLine(e, file, line, func, args...) \
    (void)(printf("ERROR: file %s, line %d, function %s:\n\tfailed assertion: (%s)\n\t", file, line, func,e), printf(args), printf("\n"), abort())

// comment out the next line if __builtin_expect causes problems
#define USE_GCC_OPTIMIZATION
#ifndef USE_GCC_OPTIMIZATION
// override system assert.h
#define dassert_(e)  \
    ((void) ((e) ? (void)0 : dassert_fileLine(#e, __FILE__, __LINE__, __func__)))
#define dassert_printf_(e, args...)  \
    ((void) ((e) ? (void)0 : dassert_printf_fileLine(#e, __FILE__, __LINE__, __func__,##args)))
#else // ifdef USE_GCC_OPTIMIZATION
// optimized version of preceding
#define dassert_(e)  \
    (__builtin_expect(!(e), 0) ? dassert_fileLine (#e, __FILE__, __LINE__, __func__) : (void)0)
#define dassert_printf(e, args...)  \
    (__builtin_expect(!(e), 0) ? dassert_printf_fileLine (#e, __FILE__, __LINE__, __func__,##args) : (void)0)
#endif // USE_GCC_OPTIMIZATION

#if(MAX_ASSERT_LEVEL>=1)
  #define assert1 dassert_
  #define assert_printf1 dassert_printf
#else
  #define assert1(e)
  #define assert_printf1(args...)
#endif
#if(MAX_ASSERT_LEVEL>=2)
  #define assert2 dassert_
  #define assert  dassert_
  #define assert_printf2 dassert_printf
  #define assert_printf  dassert_printf
#else
  #define assert2(e)
  #define assert(e)
  #define assert_printf2(args...)
  #define assert_printf(args...)
#endif
#if(MAX_ASSERT_LEVEL>=3)
  #define assert3 dassert_
  #define assert_printf3 dassert_printf
#else
  #define assert3(e)
  #define assert_printf3(args...)
#endif

// asserting specific relationships

//void assert_error_double(const char* file, int line, const char* func,
//  const char* op, const char* lhs_str, const char* rhs_str,
//  double lhs, double rhs);
#define declare_assert_errmsg(t1,t2) \
  void assert_error(const char* file, int line, const char* func, \
    const char* op, const char* lhs_str, const char* rhs_str, \
    t1 lhs, t2 rhs); 
declare_assert_errmsg(double,double); // this seems enough for all numbers
declare_assert_errmsg(int,int); // but maybe this is more efficient
declare_assert_errmsg(const char*,const char*);
// put in assert_string.h:
//#include "assert.h"
//#include<string>
//declare_assert_errmsg(const string&,const string&);

extern "C" {
int fcmp(double x1, double x2, double epsilon);
}
#ifndef USE_GCC_OPTIMIZATION
#define builtin_expect(a,b) (a)
#else
#define builtin_expect(a,b) __builtin_expect(a,b)
#endif
#define assert_not_almost_eq(lhs,rhs) \
  (fcmp(lhs, rhs, 1e-14) \
    ? (void)0 \
    : assert_error(__FILE__, __LINE__, __func__, " !=~= ", #lhs, #rhs, lhs, rhs))
#define assert_almost_eq(lhs,rhs) \
  (builtin_expect(fcmp(lhs, rhs, 1e-10),0) \
    ? assert_error(__FILE__, __LINE__, __func__, " =~= ", #lhs, #rhs, lhs, rhs) \
    : (void)0)
#define assert_divides(lhs,rhs) \
  (builtin_expect(rhs%lhs,0) \
    ? assert_error(__FILE__, __LINE__, __func__, "(divides)", #lhs, #rhs, lhs, rhs) \
    : (void)0)
#define assert_streq(lhs,rhs) \
  (builtin_expect(strcmp(lhs,rhs),0) \
    ? assert_error(__FILE__, __LINE__, __func__, "==", #lhs, #rhs, lhs, rhs) \
    : (void)0)
#define assert_op(op,lhs,rhs) \
  (builtin_expect(!(lhs op rhs),0) \
    ? assert_error(__FILE__, __LINE__, __func__, #op, #lhs, #rhs, lhs, rhs) \
    : (void)0)

// these implementations are much faster than if using 
// std::isnan(a) and std::isinf(a)
// #include <float.h> // need for DBL_MAX if invoking assert_isfinite
// check that the number is between -inf and inf
// (so that this will work even in case -ffast-math is set)
#define assert_isfinite(a) \
  (builtin_expect(!(a>=-DBL_MAX && a<=DBL_MAX),0) ? \
  (void)(printf("ERROR: file %s, line %d, function %s:\n\t: %s = %24.16e " \
    "is not finite\n", __FILE__, __LINE__, __func__,#a,a),abort()) : (void) 0)
//
// if -funsafe-math-optimizations is set this will fail to detect
// nan for old gcc compilers. We could work around this by
// implementing assert_isnum in an object module compiled without
// -funsafe-math-optimizations, but this would incur unacceptable
// function call overhead. So the user will only get checks
// against nan for modules not compiled with this option.
//
#define assert_isnum(a) \
  (builtin_expect(!(a==a),0) ? \
  (void)(printf("ERROR: file %s, line %d, function %s:\n\t: %s = %24.16e " \
    "is nan\n", __FILE__, __LINE__, __func__,#a,a),abort()) : (void) 0)

#define assert_eq(a,b) assert_op(==,a,b);
#define assert_ne(a,b) assert_op(!=,a,b);
#define assert_gt(a,b) assert_op(>,a,b);
#define assert_lt(a,b) assert_op(<,a,b);
#define assert_ge(a,b) assert_op(>=,a,b);
#define assert_le(a,b) assert_op(<=,a,b);

#endif // NDEBUG
#include <stdio.h>
#include <stdlib.h> /* for exit and abort */
#include <stdarg.h>
#include <unistd.h> /* for isatty */
#include "debug.h"

// implemented here because the new
// standard basename() in libgen.h eliminates the const
// (because standard basename returns non-const? -- seems
// like an attempted correction in the wrong direction...).
const char *dog_basename (const char *name)
{
  const char *base;
  for (base = name; *name; name++)
    if (*name == '/') base = name + 1;
  return base;
}

// === documentation of debug code ===
//
// For an example, see the test code:
//   finess/lib/tests/test_debug.cpp
//
// eprintf passes its arguments to eprintf_fileLine;
// derr causes eprintf_fileLine to be called when it goes
//   out of scope.
// 
// I have errmsg_printf_fileLine defined to call a script
//   eprintf_script which (by uncommenting a line)
//   can be configured to email the error message to the user.
// 
// To use debug, you should
//   #include "debug.h"
// which defines eprintf and dprintf.
// 
// If you insist on using streams (which I would discourage,
// since including the stream libraries increases the time
// necessary to compile an object module by over a factor of 6),
// you should
//   #include "DebugStream.h"
// to get dout (which works like cout) and
//   #include "ErrorStream.h"
// to get derr (which works like cerr, except that
// when derr goes out of scope eprintf_fileLine gets called).
// As a convenience, the line
//   #include "DebugStreams.h"
// will include both DebugStream.h and ErrorStream.h.
// 
// I have implemented three levels of debug.
// 
//   dprintf1 (or dout1) can be used for important debug,
//   dprintf2 (or dout2) can be used for standard debug, and
//   dprintf3 (or dout3) is for verbose debug
//     (such as from the inside of a loop).
//
// By default dprintf is aliased to dprintf2 and dout is
// aliased to dout2. The expectation is that these symbols are
// contextually defined per file with these default values.
// Adding the line
//   #define dout dout3
// *after* including debug.h would override this default
// 
// In debug.h notice the line
//   #define MAX_DEBUG_LEVEL 2
// If you wish to see level-3 debug you should change this to
//   #define MAX_DEBUG_LEVEL 3
// If you make this change in this file it will cause level-3
// debug to be shown universally.
// If you include the line
//   #define MAX_DEBUG_LEVEL 3
// in a particular file (*before* including debug.h)
// it will redefine the maximum debug level for that particular file.
// 
// The purpose of this constant is to provide a mechanism that
// (1) allows for file-level debug level settings and
// (2) allows debug to be removed at compile time.
// 
// The setting
//   #define MAX_DEBUG_LEVEL 2
// will cause any appearances of dprintf3 to be removed
// at compile time.  Appearances of dout3 will not be
// removed, but will be replaced with "if(0) ...",
// which (if the optimizer has any intelligence  at all)
// will incur no run-time cost.
// 
// To change the debug level universally
// at *run* time, you should call, e.g.,
// 
//   DebugLevel::set(1)
// 
// if you wish to change the debug level from its default
// value of 2.
//
// === end of documentation of debug code ===

/* To do:
   - modify to append to a log file
   - create personalized assert which email the user.
 */

int DebugLevel::level=2;  // set default debug level.

// int debug_level=2; /* default debug level */
void dprintf_fileLine(int dlevel,
  const char *func, const char *file, int line_number,
  const char *format, ...)
{
  // This test is unnecessary since the macro 
  // already takes care of it.
  //if(DebugLevel::get() < dlevel) return;
  va_list args;
  va_start(args, format);
  fprintf(stderr, "DEBUG(%d) %s(), %s:%d: ",
    dlevel, func, dog_basename(file), line_number);
  /* print out remainder of message */
  vfprintf(stderr, format, args);
  va_end(args);
  printf("\n");
}

void Wprintf_fileLine(const char *func, const char *file, int line_number,
  const char *format, ...)
{
  //if(DebugLevel::get() < dlevel) return;
  va_list args;
  va_start(args, format);
  fprintf(stderr, "WARNING, %s(), %s:%d:\n\t",
    func, dog_basename(file), line_number);
  /* print out remainder of message */
  vfprintf(stderr, format, args);
  va_end(args);
  printf("\n");
}

void eprintf_fileLine(const char *func, const char *file, int line_number,
  const char *format, ...)
{
  va_list args;
  va_start(args, format);
  fprintf(stderr, "ERROR, %s(), %s:%d:\n\t",
    func, file, line_number);
  /* print out remainder of message */
  vfprintf(stderr, format, args);
  va_end(args);
  printf("\n");
  abort();
}

/*  An example of a typical escaped string:
 *  echo "\\\"\`'()[]{}<>~@#$|\&*"
 */
static void escape_special_shell_characters(
  const char* message, char* escaped_msg)
{
  char c;
  const char* message_ptr;
  char* escaped_msg_ptr;

  message_ptr = message;
  escaped_msg_ptr = escaped_msg;

  *escaped_msg_ptr++ = '"';
  while((c=*message_ptr++) != '\0')
  {
    switch(c)
    {
      /* escape characters that are special even inside double quotes */
      case '$':
      case '`':
      case '"':
      case '\\':
        *escaped_msg_ptr++ = '\\';
    }
    *escaped_msg_ptr++ = c;
  }
  *escaped_msg_ptr++ = '"';
  *escaped_msg_ptr++ = '\0';
}


/*
 *  quit and optionally send email message
 */
void errmsg_printf_fileLine(const char *func, const char *file, int line_number,
  const char *format, ...)
{
  va_list args;
  va_start(args, format);
  fprintf(stderr, "ERROR in function %s, file %s, line %d: \n\t",
    func, file, line_number);
  /* print out remainder of message */
  vfprintf(stderr, format, args);
  va_end(args);
  // append terminating newline so user does not have to do it
  fprintf(stderr, "\n");

  // create and send message
  const int BUFFER_SIZE = 1024;
  // if the output is detached from a terminal email
  // the error message to the user.
  if(!isatty(fileno(stdout)))
  {
    char message[BUFFER_SIZE];
    char escaped_msg[2*BUFFER_SIZE+2];
    char shell_command[128+2*BUFFER_SIZE];
    char *str_pos;
    int pos;
    pos = snprintf(message, BUFFER_SIZE,
      "ERROR in function %s, file %s, line %d: \n\t",func, file, line_number);
    va_start(args, format);
    pos+=vsnprintf(message+pos,BUFFER_SIZE-pos, format, args);
    va_end(args);
    // append terminating newline to message
    snprintf(message+pos, BUFFER_SIZE-pos, "\n");

    escape_special_shell_characters(message,escaped_msg);
    str_pos = shell_command;
    str_pos+= sprintf(str_pos, 
      "${FINESS}/scripts/errmsg_printf_script ");
    sprintf(str_pos, "%s", escaped_msg);
    
    fprintf(stderr, "executing shell command: %s", shell_command);
    system(shell_command);
  }

  abort();
}

//#include "VoidStream.h"
// is this "instantiation" of this vacuous class necessary?
// I imagine not, since all of its (nonexistent) data is static
//VoidStream voidStream;

#include <iostream>
using namespace std;
#define implement_invalid_value_error(t1) \
  void invalid_value_error_fileLine(const char* file, int line, const char* func, \
    const char* type, const char* expr, t1 val) \
  { \
    std::cerr<< "ERROR in file " << file << ", line " << line  \
      << ", function " << func  \
      <<"\n\t" << type << " value: " << expr << " = " << val << endl; \
      abort(); \
  }

implement_invalid_value_error(double);
implement_invalid_value_error(int);
implement_invalid_value_error(const char*);
#include<string>
implement_invalid_value_error(const string&);

#define implement_dprintvar_fileLine(code,type) \
  void dprintvar_fileLine(int dlevel, const char* func, const char* file, int line, \
    const char* name, type val) \
  { \
    dprintf_fileLine(dlevel,func,file,line, \
      "%s == " code,name,val); \
  }
implement_dprintvar_fileLine("%s",const char*);
implement_dprintvar_fileLine("%d",int);
implement_dprintvar_fileLine("%24.16e",double);

void printvar_name_val_nl(const char*name,int val)
{ printf("%24d =: %s\n",val,name); }
void printvar_name_val_nl(const char*name,double val)
{ printf("%24.16e =: %s\n",val,name); }

void printvar_name_val(const char*name,int val)
{ printf("%s = %d",name,val); }
void printvar_name_val(const char*name,double val)
{ printf("%s = %24.16e",name,val); }

// void printvar_name_val_nl(const char*name,int val)
// { printf("%24s := %d\n",name,val); }
// void printvar_name_val_nl(const char*name,double val)
// { printf("%24s := %24.16e\n",name,val); }
#ifndef DEBUG_H
#define DEBUG_H

// on GNU systems dprintf is defined in stdio.h;
// the following line insures that this file
// is included *before* we override this symbol.
#include <stdio.h>
//////////////////////////////////////////
// for user documentation see debug.cpp //
//////////////////////////////////////////

// stuff that applies for both C-style and C++-style debug

// debug higher than this level is compiled out
// (to avoid slowing execution)
// uncomment next line to make this variable global
// #undef MAX_DEBUG_LEVEL
#ifndef MAX_DEBUG_LEVEL
  #define MAX_DEBUG_LEVEL 2
#endif
// a value of 0 turns debug completely off
#if(MAX_DEBUG_LEVEL>=0)
  #define DEBUG_ON_
#else
  #undef DEBUG_ON_
#endif

// singleton
class DebugLevel
{
private:
  static int level;
private:
  //DebugLevel(){}; // private constructor
  DebugLevel(const DebugLevel&); // prevent copy-construction
  DebugLevel& operator=(const DebugLevel&); // prevent assignment
public:
  DebugLevel(int level_){ level=level_; } // "constructor" sets (default) level.
  static void set(int level_) { level=level_; }
  static int get() { return level; }
};

// stuff specific to C-style debug

void errmsg_printf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
void eprintf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
void Wprintf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
void dprintf_fileLine(int dlevel,
                      const char *func, const char *file, int line_number, const char *format, ...);

// void eprintf_fileLine(const char* file, int line, ...);
// #define eprintf(args...) error_fileLine(__func__,__FILE__,__LINE__ , ## args)
#define errmsg_printf(args...) \
  errmsg_printf_fileLine(__func__, __FILE__, __LINE__, ## args);
#define eprintf(args...) \
  errmsg_printf_fileLine(__func__, __FILE__, __LINE__, ## args);
#define Wprintf(args...) \
  Wprintf_fileLine(__func__, __FILE__, __LINE__, ## args);
#define declare_invalid_value_error(t1) \
  void invalid_value_error_fileLine(const char* file, int line, const char* func, \
    const char* type, const char* expr, t1 val);
declare_invalid_value_error(double);
declare_invalid_value_error(int);
declare_invalid_value_error(const char*);
//#include<string>
//declare_invalid_value_error(const string&);
#define unsupported_value_error(val) invalid_value_error_fileLine( \
  __FILE__, __LINE__, __func__, "unsupported", #val, val);
#define invalid_value_error(val) invalid_value_error_fileLine( \
  __FILE__, __LINE__, __func__, "invalid", #val, val);

const char *dog_basename (const char *name);

inline bool debug_printFileLine(int dlevel, const char* func,
  const char *file, int line_number)
{
  if(DebugLevel::get() >= dlevel)
  {
    fprintf(stderr, "DEBUG(%d) %s(), %s:%d: ",
      dlevel, func, dog_basename(file), line_number);
    return true;
  }
  return false;
}

// version of dprintfn that returns true

// compile out unused debug
#ifndef DEBUG_ON
#define DEBUG_ON 1
#endif
#if(DEBUG_ON!=0)
  #define dprintfn(dlevel, args...) if(DebugLevel::get() >= dlevel) \
    dprintf_fileLine(dlevel,__func__,  __FILE__, __LINE__, ## args);
  #define debugn(dlevel) debug_printFileLine(dlevel, __func__, __FILE__, __LINE__)
  #define dbprintfn(dlevel, str, args...) debugn(dlevel) && \
    printf(str, ## args), printf("\n"), true
#else
  #define dprintfn(args...) 
  #define dbprintfn(args...) false
  #define debugn(dlevel) false
#endif
#if(MAX_DEBUG_LEVEL>=1)
  #define dprintf1(args...) if(DebugLevel::get() >= 1) \
    dprintf_fileLine(1, __func__, __FILE__, __LINE__, ## args);
  #define dprint1(var) if(DebugLevel::get() >= 1) \
    dprintvar_fileLine(1,__func__,__FILE__,__LINE__,#var,var);
  #define debug1 debug_printFileLine(1, __func__, __FILE__, __LINE__)
#else
  #define dprintf1(args...) 
  #define dprint1(args...) 
  #define debug1 false
#endif
#if(MAX_DEBUG_LEVEL>=2)
  #define dprintf2(args...) if(DebugLevel::get() >= 2) \
    dprintf_fileLine(2, __func__, __FILE__, __LINE__, ## args);
  #define dprint2(var) if(DebugLevel::get() >= 2) \
    dprintvar_fileLine(2,__func__,__FILE__,__LINE__,#var,var);
  #define debug2 debug_printFileLine(2, __func__, __FILE__, __LINE__)
#else
  #define dprintf2(args...) 
  #define dprint2(args...) 
  #define debug2 false
#endif
#if(MAX_DEBUG_LEVEL>=3)
  #define dprintf3(args...) if(DebugLevel::get() >= 3) \
    dprintf_fileLine(3, __func__, __FILE__, __LINE__, ## args);
  #define dprint3(var) if(DebugLevel::get() >= 3) \
    dprintvar_fileLine(3,__func__,__FILE__,__LINE__,#var,var);
  #define debug3 debug_printFileLine(3, __func__, __FILE__, __LINE__)
#else
  #define dprintf3(args...) 
  #define dprint3(args...) 
  #define debug3 false
#endif

// debug level of dprintf and dprint is intended to be contextually defined;
// define a default debug level for dprintf
// (user is expected to override this definition)
// 
#define dprintf dprintf2
#define dprint dprint2

#define declare_dprintvar_fileLine(type) \
  void dprintvar_fileLine(int,const char*,const char*,int,const char*,type);
declare_dprintvar_fileLine(int);
declare_dprintvar_fileLine(double);
declare_dprintvar_fileLine(const char*);

#define println(arg, args...) \
  printf(arg, ## args); printf("\n");
void printvar_name_val(const char*name,int val);
void printvar_name_val(const char*name,double val);
void printvar_name_val_nl(const char*name,int val);
void printvar_name_val_nl(const char*name,double val);
#define printvn(var) printvar_name_val_nl(#var,var);
#define printv(var) printvar_name_val(#var,var);
#define printvc(var) printvar_name_val(#var,var); printf(", ");

#endif
#ifndef _UTIL_H_
#define _UTIL_H_

#include <sstream>
template<typename T>
T stringToAny(const std::string& s){
    T ret;
    std::istringstream iss(s);
    iss >> std::boolalpha >> ret;
    return ret;
}

template<typename T>
std::string anyToString(const T& a){
    std::stringstream ss;
    ss << std::boolalpha << a;
    return ss.str();
}

#include <cstdio>
inline bool existFile(const std::string& filename){
    using std::FILE;
    using std::fopen;
    using std::fclose;
    if(FILE *file = fopen(filename.c_str(), "r")){
	fclose(file);
	return true;
    } else{
	return false;
    }
}

#include <iostream>
#include <cstdlib>
#include <stdexcept>
//#define TERMINATE_ABORT
#undef TERMINATE_ABORT
inline void terminate(const std::string& message){
#ifdef TERMINATE_ABORT
    std::cerr << message << std::endl;
    std::terminate();
#else
    throw std::runtime_error(message);
#endif
}

inline std::string& remove_trailing_spaces(std::string& s){
    std::size_t found = s.find_last_not_of(" \t");
    if(found != std::string::npos)
        s.erase(found + 1);
    else
        s.clear();
    return s;
}

inline std::string read_entire_stream(std::istream& is){
    std::string str((std::istreambuf_iterator<char>(is)),
               std::istreambuf_iterator<char>());
    return str;
}

#endif
// turn on CHECK_BOUNDS before including tensors.h
// so that the non-inline versions of get and set are declared.
// This way they are implemented in this file so that
// elsewhere in the code one can choose whether to use
// the in-line or bounds-checking versions of get and set
#ifndef CHECK_BOUNDS
#define CHECK_BOUNDS
#endif
#include "tensors.h"
#include "debug.h"
#include "assert.h"
#include<new>
#define bounds_check() expensive_bounds_check();
//#define bounds_check() cheap_bounds_check();
#define PARALLEL_SIZE 512
#include <limits> // for numeric_limits

// class iTensor2

iTensor2::iTensor2(int numRows, int numColumns)
{
    assert_ge(numRows,0);
    assert_ge(numColumns,0);

    size = numRows*numColumns;
    vec = new int[size];
    rows = numRows;
    columns = numColumns;
}

iTensor2::iTensor2(const iTensor2& anotheriTensor2)
{
    rows = anotheriTensor2.rows;
    columns = anotheriTensor2.columns;
    size = anotheriTensor2.size;
    vec = new int[size];

    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++)
        { vec[i] = anotheriTensor2.vec[i]; }
    else
#pragma omp parallel for
        for (int i=0; i<size; i++)
        { vec[i] = anotheriTensor2.vec[i]; }
}

iTensor2::~iTensor2()
{
    delete[] vec;
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_gt(n1,0); assert_le(n1,rows); \
assert_gt(n2,0); assert_le(n2,columns); 
#define cheap_bounds_check() \
    assert_printf(k<size && k>=0, "n1=%d, n2=%d",n1,n2);

const int& iTensor2::get(int n1, int n2) const
{
    int k = (n1-1)*columns + (n2-1);
    bounds_check();
    return vec[k];
}

void iTensor2::set(int n1, int n2, int value)
{
    int k = (n1-1)*columns + (n2-1);
    bounds_check();
    vec[k] = value;
}

// class iTensor3

iTensor3::iTensor3(int n1, int n2, int n3)
{
    assert_printf(n1>0 && n2>0 && n3>0, "n1=%d, n2=%d, n3=%d",n1,n2,n3);
    size = n1*n2*n3;
    vec = new int[size];
    numElem1 = n1;
    numElem2 = n2;
    numElem3 = n3;
}

iTensor3::iTensor3(const iTensor3& anotheriTensor3)
{
    numElem1 = anotheriTensor3.numElem1;
    numElem2 = anotheriTensor3.numElem2;
    numElem3 = anotheriTensor3.numElem3;
    size = anotheriTensor3.size;
    vec = new int[size];

    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++)
        { vec[i] = anotheriTensor3.vec[i]; }
    else
#pragma omp parallel for
        for (int i=0; i<size; i++)
        { vec[i] = anotheriTensor3.vec[i]; }
}

iTensor3::~iTensor3()
{
    delete[] vec;
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_gt(n1,0); assert_le(n1,numElem1); \
assert_gt(n2,0); assert_le(n2,numElem2); \
assert_gt(n3,0); assert_le(n3,numElem3); 
#define cheap_bounds_check() \
    assert_printf(k<size && k>=0, "n1=%d, n2=%d, n3=%d",n1,n2,n3);

const int& iTensor3::get(int n1, int n2, int n3) const
{
    int k = ((n1-1)*numElem2 + (n2-1))*numElem3 + (n3-1);
    bounds_check();
    return vec[k];
}

void iTensor3::set(int n1, int n2, int n3, int value)
{
    int k = ((n1-1)*numElem2 + (n2-1))*numElem3 + (n3-1);
    bounds_check();
    vec[k] = value;
}

// === major section: itensors (arrays of int) ===

// class iTensorBase

void iTensorBase::init()
{
    assert_printf(size>=0, "size=%d", size);
    vec = new int[size];
}

iTensorBase::iTensorBase(const iTensorBase& in) :
    size(in.size)
{
    vec = new int[size];
    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++) vec[i] = in.vec[i];
    else
#pragma omp parallel for
        for (int i=0; i<size; i++) vec[i] = in.vec[i];
}

iTensorBase::~iTensorBase()
{
    delete[] vec;
}

const int& iTensorBase::vget(int k) const
{
    assert_ge(k,0); assert_lt(k,size);
    return vec[k];
}

int& iTensorBase::vfetch(int k)
{
    assert_ge(k,0); assert_lt(k,size);
    return vec[k];
}

void iTensorBase::vset(int k, int value)
{
    assert_ge(k,0); assert_lt(k,size);
    vec[k]=value;
}

void iTensorBase::setall(int value)
{
    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++) vec[i] = value;
    else
#pragma omp parallel for
        for (int i=0; i<size; i++) vec[i] = value;
}

// === major section: dtensors (arrays of double) ===

// class dTensorBase

void dTensorBase::init()
{
    assert_printf(size>=0, "size=%d", size);
    vec = new double[size];
#ifdef CHECK_INIT
    setall(std::numeric_limits<double>::quiet_NaN());
    // setall(0./0); // initialize all entries to nan
#endif
}

dTensorBase::dTensorBase(const dTensorBase& in, CopyMode::Enum copyMode) :
    size(in.size)
{
    switch(copyMode)
    {
        // should make vec counted_ptr first
        //case CopyMode::SHALLOW:
        //  vec(in.vec);
        //  break;
        case CopyMode::DIMS:
            vec = new double[size];
            break;
        case CopyMode::DEEP:
            vec = new double[size];
            copyfrom(in);
            break;
        default:
            unsupported_value_error(copyMode);
    }
}

void dTensorBase::copyfrom(const dTensorBase& in)
{
    if(this!=&in)
    {
        assert_eq(size,in.size);
        if(size<PARALLEL_SIZE)
            for (int i=0; i<size; i++) vec[i] = in.vec[i];
        else
#pragma omp parallel for
            for (int i=0; i<size; i++) vec[i] = in.vec[i];
    }
}

dTensorBase::~dTensorBase()
{
    delete[] vec;
}

const double& dTensorBase::vget(int k) const
{
    assert_ge(k,0); assert_lt(k,size);
#ifdef CHECK_INIT
    if(!(vec[k]==vec[k])) eprintf("vec[%d]=%24.16e",k,vec[k]);
    // double ret = vec[k];
    // assert_eq(ret,ret);
#endif

    return vec[k];
}

double& dTensorBase::vfetch(int k)
{
    assert_ge(k,0); assert_lt(k,size);
#ifdef CHECK_INIT
    if(!(vec[k]==vec[k])) eprintf("vec[%d]=%24.16e",k,vec[k]);
#endif
    return vec[k];
}

void dTensorBase::vset(int k, double value)
{
    assert_ge(k,0); assert_lt(k,size);
#ifdef CHECK_INIT
    if(!(value==value)) eprintf("for k=%d value=%24.16e",k,value);
#endif
    vec[k]=value;
}

void dTensorBase::setall(double value)
{
    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++) vec[i] = value;
    else
#pragma omp parallel for
        for (int i=0; i<size; i++) vec[i] = value;
}

bool dTensorBase::check()
{
    if(size<PARALLEL_SIZE)
        for (int i=0; i<size; i++){
            assert_printf(vec[i]==vec[i],"vec[%d]=%24.16e",i,vec[i]);
        }
    else
#pragma omp parallel for
        for (int i=0; i<size; i++){
            assert_printf(vec[i]==vec[i],"vec[%d]=%24.16e",i,vec[i]);
        }
}

// --- section: multidimensional tensor base classes ---

// class dTensor5d

// check parameters and update derived parameters.
void dTensor5d::init()
{
    assert_printf(s1>=0 && s2>=0 && s3>=0 && s4>=0 && s5>=0,
            "b1=%d, b2=%d, b3=%d, b4=%d, b5=%d"
            "s1=%d, s2=%d, s3=%d, s4=%d, s5=%d",
            b1,b2,b3,b4,b5,
            s1,s2,s3,s4,s5);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    e3 = b3 + s3 - 1;
    e4 = b4 + s4 - 1;
    e5 = b5 + s5 - 1;
    size = s1*s2*s3*s4*s5;
    dTensorBase::init();
}

dTensor5d::dTensor5d(
        int s1i, int s2i, int s3i, int s4i, int s5i,
        int b1i, int b2i, int b3i, int b4i, int b5i ) :
    s1(s1i), s2(s2i), s3(s3i), s4(s4i), s5(s5i), 
    b1(b1i), b2(b2i), b3(b3i), b4(b4i), b5(b5i)
{ init(); }

dTensor5d::dTensor5d(const dTensor5d& in, CopyMode::Enum copyMode) 
    : dTensorBase(in, copyMode),
    b1(in.b1), b2(in.b2), b3(in.b3), b4(in.b4), b5(in.b5),
    e1(in.e1), e2(in.e2), e3(in.e3), e4(in.e4), e5(in.e5),
    s1(in.s1), s2(in.s2), s3(in.s3), s4(in.s4), s5(in.s5) { }

    // alias for assert_eq
    //#define ae(a,b) assert_op(==,a,b);
#define ae assert_eq
void dTensor5d::copyfrom(const dTensor5d& in)
{
    ae(b1,in.b1); ae(b2,in.b2); ae(b3,in.b3); ae(b4,in.b4); ae(b5,in.b5);
    ae(e1,in.e1); ae(e2,in.e2); ae(e3,in.e3); ae(e4,in.e4); ae(e5,in.e5);
    ae(s1,in.s1); ae(s2,in.s2); ae(s3,in.s3); ae(s4,in.s4); ae(s5,in.s5);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2); \
assert_ge(n3,b3); assert_le(n3,e3); \
assert_ge(n4,b4); assert_le(n4,e4); \
assert_ge(n5,b5); assert_le(n5,e5); 
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d, n3=%d, n4=%d, n5=%d", n1, n2, n3, n4, n5);

const double& dTensor5d::get(int n1,int n2,int n3,int n4,int n5) const
{
    int k = getidx(n1,n2,n3,n4,n5);
    bounds_check();
    return vec[k];
}

double& dTensor5d::fetch(int n1,int n2,int n3,int n4,int n5)
{
    int k = getidx(n1,n2,n3,n4,n5);
    bounds_check();
    return vec[k];
}


void dTensor5d::set(int n1,int n2,int n3,int n4,int n5,double value)
{
    int k = getidx(n1,n2,n3,n4,n5);
    bounds_check();
    vec[k] = value;
}

// class dTensor4d

// check parameters and update derived parameters.
void dTensor4d::init()
{
    assert_printf(s1>=0 && s2>=0 && s3>=0 && s4>=0,
            "b1=%d, b2=%d, b3=%d, b4=%d"
            "s1=%d, s2=%d, s3=%d, s4=%d",
            b1,b2,b3,b4,
            s1,s2,s3,s4);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    e3 = b3 + s3 - 1;
    e4 = b4 + s4 - 1;
    size = s1*s2*s3*s4;
    dTensorBase::init();
}

dTensor4d::dTensor4d(
        int s1i, int s2i, int s3i, int s4i,
        int b1i, int b2i, int b3i, int b4i ) :
    s1(s1i), s2(s2i), s3(s3i), s4(s4i), 
    b1(b1i), b2(b2i), b3(b3i), b4(b4i)
{ init(); }

dTensor4d::dTensor4d(const dTensor4d& in, CopyMode::Enum copyMode)
    : dTensorBase(in, copyMode),
    b1(in.b1), b2(in.b2), b3(in.b3), b4(in.b4),
    e1(in.e1), e2(in.e2), e3(in.e3), e4(in.e4),
    s1(in.s1), s2(in.s2), s3(in.s3), s4(in.s4) { }

void dTensor4d::copyfrom(const dTensor4d& in)
{
    ae(b1,in.b1); ae(b2,in.b2); ae(b3,in.b3); ae(b4,in.b4);
    ae(e1,in.e1); ae(e2,in.e2); ae(e3,in.e3); ae(e4,in.e4);
    ae(s1,in.s1); ae(s2,in.s2); ae(s3,in.s3); ae(s4,in.s4);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2); \
assert_ge(n3,b3); assert_le(n3,e3); \
assert_ge(n4,b4); assert_le(n4,e4);
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d, n3=%d, n4=%d", n1, n2, n3, n4);

const double& dTensor4d::get(int n1,int n2,int n3,int n4) const
{
    int k = getidx(n1,n2,n3,n4);
    bounds_check();
    return vec[k];
}

double& dTensor4d::fetch(int n1,int n2,int n3,int n4)
{
    int k = getidx(n1,n2,n3,n4);
    bounds_check();
    return vec[k];
}

void dTensor4d::set(int n1,int n2,int n3,int n4,double value)
{
    int k = getidx(n1,n2,n3,n4);
    bounds_check();
    vec[k] = value;
}

// class dTensor3d

// check parameters and update derived parameters.
void dTensor3d::init()
{
    assert_printf(s1>=0 && s2>=0 && s3>=0,
            "b1=%d, b2=%d, b3=%d"
            "s1=%d, s2=%d, s3=%d",
            b1,b2,b3,
            s1,s2,s3);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    e3 = b3 + s3 - 1;
    size = s1*s2*s3;
    dTensorBase::init();
}

dTensor3d::dTensor3d(
        int s1i, int s2i, int s3i,
        int b1i, int b2i, int b3i ) :
    s1(s1i), s2(s2i), s3(s3i), 
    b1(b1i), b2(b2i), b3(b3i)
{ init(); }

dTensor3d::dTensor3d(const dTensor3d& in) : dTensorBase(in),
    b1(in.b1), b2(in.b2), b3(in.b3),
    e1(in.e1), e2(in.e2), e3(in.e3),
    s1(in.s1), s2(in.s2), s3(in.s3) { }

void dTensor3d::copyfrom(const dTensor3d& in)
{
    ae(b1,in.b1); ae(b2,in.b2); ae(b3,in.b3);
    ae(e1,in.e1); ae(e2,in.e2); ae(e3,in.e3);
    ae(s1,in.s1); ae(s2,in.s2); ae(s3,in.s3);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2); \
assert_ge(n3,b3); assert_le(n3,e3);
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d, n3=%d", n1, n2, n3);

const double& dTensor3d::get(int n1,int n2,int n3) const
{
    int k = getidx(n1,n2,n3);
    bounds_check();
    return vec[k];
}

void dTensor3d::set(int n1,int n2,int n3,double value)
{
    int k = getidx(n1,n2,n3);
    bounds_check();
    vec[k] = value;
}

// class dTensor2d

// check parameters and update derived parameters.
void dTensor2d::init()
{
    assert_printf(s1>=0 && s2>=0,
            "b1=%d, b2=%d"
            "s1=%d, s2=%d",
            b1,b2,
            s1,s2);
    e1 = b1 + s1 - 1;
    e2 = b2 + s2 - 1;
    size = s1*s2;
    dTensorBase::init();
}

dTensor2d::dTensor2d(
        int s1i, int s2i,
        int b1i, int b2i ) :
    s1(s1i), s2(s2i), 
    b1(b1i), b2(b2i)
{ init(); }

dTensor2d::dTensor2d(const dTensor2d& in) : dTensorBase(in),
    b1(in.b1), b2(in.b2),
    e1(in.e1), e2(in.e2),
    s1(in.s1), s2(in.s2) { }

void dTensor2d::copyfrom(const dTensor2d& in)
{
    ae(b1,in.b1); ae(b2,in.b2);
    ae(e1,in.e1); ae(e2,in.e2);
    ae(s1,in.s1); ae(s2,in.s2);
    dTensorBase::copyfrom(in);
}

#undef expensive_bounds_check
#undef cheap_bounds_check
#define expensive_bounds_check() \
    assert_ge(n1,b1); assert_le(n1,e1); \
assert_ge(n2,b2); assert_le(n2,e2);
#define cheap_bounds_check() assert_printf(k<size && k>=0, \
        "n1=%d, n2=%d", n1, n2);

const double& dTensor2d::get(int n1,int n2) const
{
    int k = getidx(n1,n2);
    bounds_check();
    return vec[k];
}

double& dTensor2d::fetch(int n1,int n2)
{
    int k = getidx(n1,n2);
    bounds_check();
    return vec[k];
}

void dTensor2d::set(int n1,int n2,double value)
{
    int k = getidx(n1,n2);
    bounds_check();
    vec[k] = value;
}

void dTensor1d::copyfrom(const dTensor1d& in)
{
    ae(b1,in.b1);
    dTensorBase::copyfrom(in);
}

// --- section: multidimensional 1-based tensor classes ---

// class dTensorBC5

void dTensorBC5::copyfrom(const dTensorBC5& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor5d::copyfrom(in);
}

dTensorBC5::dTensorBC5(const dTensorBC5& in, CopyMode::Enum copyMode)
: dTensor5d(in,copyMode),
  S1(in.S1), S2(in.S2), S3(in.S3), S4(in.S4), S5(in.S5),
  mbc(in.mbc), ndims(in.ndims)
{ }

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC5::dTensorBC5(int S1i, int S2i, int S3i, int S4i, int S5i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), S3(S3i), S4(S4i), S5(S5i), 
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i; s3=S3i; s4=S4i; s5=S5i; 
    b1=1;   b2=1;   b3=1;   b4=1;   b5=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 5: b5=1-mbc; s5=S5+2*mbc;
        case 4: b4=1-mbc; s4=S4+2*mbc;
        case 3: b3=1-mbc; s3=S3+2*mbc;
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}

// class dTensorBC4

void dTensorBC4::copyfrom(const dTensorBC4& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor4d::copyfrom(in);
}

dTensorBC4::dTensorBC4(const dTensorBC4& in, CopyMode::Enum copyMode)
: dTensor4d(in,copyMode),
  S1(in.S1), S2(in.S2), S3(in.S3), S4(in.S4), 
  mbc(in.mbc), ndims(in.ndims)
{ }

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC4::dTensorBC4(int S1i, int S2i, int S3i, int S4i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), S3(S3i), S4(S4i), 
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i; s3=S3i; s4=S4i;
    b1=1;   b2=1;   b3=1;   b4=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 4: b4=1-mbc; s4=S4+2*mbc;
        case 3: b3=1-mbc; s3=S3+2*mbc;
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}

// class dTensorBC3

void dTensorBC3::copyfrom(const dTensorBC3& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor3d::copyfrom(in);
}

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC3::dTensorBC3(int S1i, int S2i, int S3i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), S3(S3i), 
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i; s3=S3i;
    b1=1;   b2=1;   b3=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 3: b3=1-mbc; s3=S3+2*mbc;
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}

// class dTensorBC2

void dTensorBC2::copyfrom(const dTensorBC2& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor2d::copyfrom(in);
}

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC2::dTensorBC2(int S1i, int S2i,
        int mbcin, int ndimsin) :
    S1(S1i), S2(S2i), 
    mbc(mbcin), ndims(ndimsin)
{
    s1=S1i; s2=S2i;
    b1=1;   b2=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 2: b2=1-mbc; s2=S2+2*mbc;
        case 1: b1=1-mbc; s1=S1+2*mbc;
        case 0: ;
    }
    init();
}

// class dTensorBC1

void dTensorBC1::copyfrom(const dTensorBC1& in)
{
    assert_eq(mbc,in.mbc); assert_eq(ndims,in.ndims);
    dTensor1d::copyfrom(in);
}

// ndims: number of dimensions that have mbc layers of boundary cells
dTensorBC1::dTensorBC1(int S1i,
        int mbcin, int ndimsin) :
    S1(S1i),
    mbc(mbcin), ndims(ndimsin)
{
    size=S1i;
    b1=1;
    assert_ge(mbc,0);
    switch(ndims)
    {
        default: invalid_value_error(ndims);
        case 1: b1=1-mbc; size=S1+2*mbc;
        case 0: ;
    }
    init();
}

#ifndef _TENSORS_H_
#define _TENSORS_H_
#include "tensors1d.h"
//#include "debug.h" // for invalid_value_error

// --------------------------------------------------------------------------
// tensors.h defines DoGPack's multidimensional array classes.
//   Array indices increment in odometer order.
//
// [id]TensorBC[1-9] are tensor classes whose first two components
//   represent mesh coordinates; mbc is the number of layers of ghost cells
//   at the boundary.
//
// For example:
//
// In dogpack code the 3d state is declared as "dTensorBC5 q" and
//   q(i,j,k,m,ell) represents the state at x-index i, y-index j, 
//   and z-index k for equation index m and for polynomial basis index ell.
//   Values of i between 1 and q.getsize(1),
//   values of j between 1 and q.getsize(2), and
//   values of k between 1 and q.getsize(3) represent state cells.
//   Values of i between 1-mbc and 0
//          and between q.getsize(1)+1 and q.getsize(1)+q.getmbc(),
//   values of j between 1-mbc and 0
//          and between q.getsize(2)+1 and q.getsize(2)+q.getmbc(), and
//   values of k between 1-mbc and 0
//          and between q.getsize(3)+1 and q.getsize(3)+q.getmbc()
//     represent ghost cell values.
//
//   mbc is 2 in most applications.  Specifically:
//   - Advancing boundary cells requires 1st layer of ghost cells.
//   - Limiting advanced values of boundary cells (which is needed
//       for stability if you want higher than first order accuracy
//       in space) requires advancing 1st layer of ghost cells.
//   - Advancing 1st layer of ghost cells requires 2nd layer of ghost cells.
// --------------------------------------------------------------------------
//
// By default get and set methods are defined
//   inline without bounds checking.
// To perform bounds checking, compile with -DCHECK_BOUNDS
//   You can put #undef CHECK_BOUNDS immediately prior to
//   #include "tensors.h" in files that have already been
//   debugged.)

class iTensor2
{
    public:
        iTensor2(int n1, int n2);
        // Constructor
        // POST: Create a matrix with n1 rows and n2 columns

        iTensor2(const iTensor2& anotheriTensor2);
        // Copy constructor
        // POST: New tensor created with size and contents same as anotheriTensor2

        ~iTensor2();
        // Destructor
        // POST: iTensor no longer exists

#ifndef CHECK_BOUNDS
        const int& get(int n1, int n2) const {
            return vec[ (n1-1)*columns + (n2-1) ]; }
        void set(int n1, int n2, int value) {
            vec[ (n1-1)*columns + (n2-1) ] = value; }
#else
        const int& get(int n1, int n2) const;
        // POST: Get (n1,n2)^(th) entry in tensor

        void set(int n1, int n2, int value);
        // POST: Set (n1,n2)^(th) entry in tensor to "value"
#endif

        int getsize(int n) const
            // POST: if n==1: returns number of rows
            //       if n==2: returns number of columns
        {
            switch(n)
            {
                case 1: return rows;
                case 2: return columns;
                default: return 1;
            }
        }

    private:
        int* vec;
        int rows,columns;
        int size;
};

class iTensor3
{
    public:
        iTensor3(int n1, int n2, int n3);
        // Constructor
        // POST: Create a tensor of size (n1,n2,n3)

        iTensor3(const iTensor3& anotheriTensor3);
        // Copy constructor
        // POST: New tensor created with size and contents same as anotheriTensor3

        ~iTensor3();
        // Destructor
        // POST: iTensor no longer exists

#ifndef CHECK_BOUNDS
        const int& get(int n1, int n2, int n3) const {
            return vec[ ((n1-1)*numElem2 + (n2-1))*numElem3 + (n3-1) ]; }
        void set(int n1, int n2, int n3, int value) {
            vec[ ((n1-1)*numElem2 + (n2-1))*numElem3 + (n3-1) ] = value; }
#else
        const int& get(int n1, int n2, int n3) const;
        // POST: Get (n1,n2,n3)^(th) entry in tensor

        void set(int n1, int n2, int n3, int value);
        // POST: Set (n1,n2,n3)^(th) entry in tensor to "value"
#endif

        int getsize(int n) const
            // POST: if n==1: returns number of elements in first index
            //       if n==2: returns number of elements in second index
            //       if n==3: returns number of elements in third index
        {
            switch(n)
            {
                case 1: return numElem1;
                case 2: return numElem2;
                case 3: return numElem3;
                default: 1;
            }
        }

    private:
        int* vec;
        int numElem1,numElem2,numElem3;
        int size;
};

// === major section: dtensors (arrays of double) ===

// --- section: multidimensional tensor base classes ---

class dTensor5d : public dTensorBase
{
    // data
    protected:
        int s1, s2, s3, s4, s5;
        int b1, b2, b3, b4, b5;
        int e1, e2, e3, e4, e5;
        // methods
    private: // disabled
        dTensor5d& operator=(const dTensor5d& in);
    protected:
        dTensor5d(){};
        void init();
    public:
        int getidx(int n1, int n2, int n3, int n4, int n5) const
        {
            int k = n1-b1;
            k *= s2; k += (n2-b2);
            k *= s3; k += (n3-b3);
            k *= s4; k += (n4-b4);
            k *= s5; k += (n5-b5);
            return k;
        }
        int getsize(int n) const
        {
            switch(n)
            {
                case 5: return s5;
                case 4: return s4;
                case 3: return s3;
                case 2: return s2;
                case 1: return s1;
                default: return 1;
            }
        }
    public:
        // constructor takes size and initial index in each dimension
        dTensor5d(
                int s1i, int s2i, int s3i, int s4i, int s5i,
                int b1i, int b2i, int b3i, int b4i, int b5i );
        dTensor5d(const dTensor5d& in, CopyMode::Enum copyMode);
        void copyfrom(const dTensor5d& in);
#ifdef CHECK_BOUNDS
        const double& get(int n1,int n2,int n3,int n4,int n5) const;
        double& fetch(int n1,int n2,int n3,int n4,int n5);
        void set(int n1,int n2,int n3,int n4,int n5, double value);
#else
        const double& get(int n1,int n2,int n3,int n4,int n5) const
        { return vec[getidx(n1,n2,n3,n4,n5)]; }
        double& fetch(int n1,int n2,int n3,int n4,int n5)
        { return vec[getidx(n1,n2,n3,n4,n5)]; }
        void set(int n1,int n2,int n3,int n4,int n5, double value)
        { vec[getidx(n1,n2,n3,n4,n5)] = value; }
#endif
};

class dTensor4d : public dTensorBase
{
    // data
    protected:
        int s1, s2, s3, s4;
        int b1, b2, b3, b4;
        int e1, e2, e3, e4;
        // methods
    private: // disabled
        dTensor4d& operator=(const dTensor4d& in);
    protected:
        dTensor4d(){};
        dTensor4d(const dTensor4d& in, CopyMode::Enum copyMode=CopyMode::DEEP);
        void init();
    public:
        int getidx(int n1, int n2, int n3, int n4) const
        {
            int k = n1-b1;
            k *= s2; k += (n2-b2);
            k *= s3; k += (n3-b3);
            k *= s4; k += (n4-b4);
            return k;
        }
        int getsize(int n) const
        {
            switch(n)
            {
                case 4: return s4;
                case 3: return s3;
                case 2: return s2;
                case 1: return s1;
                default: return 1;
            }
        }
    public:
        // constructor takes size and initial index in each dimension
        dTensor4d(
                int s1i, int s2i, int s3i, int s4i,
                int b1i, int b2i, int b3i, int b4i );
        void copyfrom(const dTensor4d& in);
#ifdef CHECK_BOUNDS
        const double& get(int n1,int n2,int n3,int n4) const;
        double& fetch(int n1,int n2,int n3,int n4);
        void set(int n1,int n2,int n3,int n4, double value);
#else
        const double& get(int n1,int n2,int n3,int n4) const
        { return vec[getidx(n1,n2,n3,n4)]; }
        double& fetch(int n1,int n2,int n3,int n4)
        { return vec[getidx(n1,n2,n3,n4)]; }
        void set(int n1,int n2,int n3,int n4, double value)
        { vec[getidx(n1,n2,n3,n4)] = value; }
#endif
};

class dTensor3d : public dTensorBase
{
    // data
    protected:
        int s1, s2, s3;
        int b1, b2, b3;
        int e1, e2, e3;
        // methods
    private: // disabled
        dTensor3d& operator=(const dTensor3d& in);
    protected:
        dTensor3d(){};
        void init();
    public:
        int getidx(int n1, int n2, int n3) const
        {
            int k = n1-b1;
            k *= s2; k += (n2-b2);
            k *= s3; k += (n3-b3);
            return k;
        }
        int getsize(int n) const
        {
            switch(n)
            {
                case 3: return s3;
                case 2: return s2;
                case 1: return s1;
                default: return 1;
            }
        }
    public:
        // constructor takes size and initial index in each dimension
        dTensor3d(
                int s1i, int s2i, int s3i,
                int b1i, int b2i, int b3i );
        dTensor3d(const dTensor3d& in);
        void copyfrom(const dTensor3d& in);
#ifdef CHECK_BOUNDS
        const double& get(int n1,int n2,int n3) const;
        void set(int n1,int n2,int n3, double value);
#else
        const double& get(int n1,int n2,int n3) const
        { return vec[getidx(n1,n2,n3)]; }
        void set(int n1,int n2,int n3, double value)
        { vec[getidx(n1,n2,n3)] = value; }
#endif
};

class dTensor2d : public dTensorBase
{
    // data
    protected:
        int s1, s2;
        int b1, b2;
        int e1, e2;
        // methods
    private: // disabled
        dTensor2d& operator=(const dTensor2d& in);
    protected:
        dTensor2d(){};
        void init();
    public:
        int getidx(int n1, int n2) const
        {
            int k = n1-b1;
            k *= s2; k += (n2-b2);
            return k;
        }
        int getsize(int n) const
        {
            switch(n)
            {
                case 2: return s2;
                case 1: return s1;
                default: return 1;
            }
        }
    public:
        // constructor takes size and initial index in each dimension
        dTensor2d(
                int s1i, int s2i,
                int b1i, int b2i );
        dTensor2d(const dTensor2d& in);
        void copyfrom(const dTensor2d& in);
#ifdef CHECK_BOUNDS
        const double& get(int n1,int n2) const;
        double& fetch(int n1,int n2);
        void set(int n1,int n2, double value);
#else
        const double& get(int n1,int n2) const
        { return vec[getidx(n1,n2)]; }
        double& fetch(int n1,int n2)
        { return vec[getidx(n1,n2)]; }
        void set(int n1,int n2, double value)
        { vec[getidx(n1,n2)] = value; }
#endif
};

// --- section: multidimensional 1-based tensor classes ---

class dTensor5 : public dTensor5d
{
    private: // disabled
        dTensor5& operator=(const dTensor5& in){ copyfrom(in); return *this; }
    public:
        dTensor5(int s1i, int s2i, int s3i, int s4i, int s5i) :
            dTensor5d(s1i,s2i,s3i,s4i,s5i,1,1,1,1,1) { }
        void copyfrom(const dTensor5& in){ dTensor5d::copyfrom(in); }

        // For speed we override the defaults
        // We can delete this; it gives no detectable speedup
#ifndef CHECK_BOUNDS
        int getidx(int n1, int n2, int n3, int n4, int n5) const
        {
            int k = n1-1;
            k *= s2; k += (n2-1);
            k *= s3; k += (n3-1);
            k *= s4; k += (n4-1);
            k *= s5; k += (n5-1);
            return k;
        }
        const double& get(int n1,int n2,int n3,int n4,int n5) const
        { return vec[getidx(n1,n2,n3,n4,n5)]; }
        void set(int n1,int n2,int n3,int n4,int n5, double value)
        { vec[getidx(n1,n2,n3,n4,n5)] = value; }
#endif
};

class dTensor4 : public dTensor4d
{
    public:
        dTensor4(int s1i, int s2i, int s3i, int s4i) :
            dTensor4d(s1i,s2i,s3i,s4i,1,1,1,1) { }
        void copyfrom(const dTensor4& in){ dTensor4d::copyfrom(in); }
    private: // disabled
        dTensor4& operator=(const dTensor4& in){ copyfrom(in); return *this; }
};

class dTensor3 : public dTensor3d
{
    public:
        dTensor3(int s1i, int s2i, int s3i) :
            dTensor3d(s1i,s2i,s3i,1,1,1) { }
        void copyfrom(const dTensor3& in){ dTensor3d::copyfrom(in); }
    private: // disabled
        dTensor3& operator=(const dTensor3& in){ copyfrom(in); return *this; }
};

class dTensor2 : public dTensor2d
{
    public:
        dTensor2(int s1i, int s2i) :
            dTensor2d(s1i,s2i,1,1) { }
        void copyfrom(const dTensor2& in){ dTensor2d::copyfrom(in); }
    private: // disabled
        dTensor2& operator=(const dTensor2& in){ copyfrom(in); return *this; }
};

// --- section: multidimensional boundary condition (BC) tensor classes ---

class dTensorBC5 : public dTensor5d
{
    private: // disabled
        dTensorBC5& operator=(const dTensorBC5& in);
    public:
        dTensorBC5(int s1i, int s2i, int s3i, int s4i, int s5i,
                int mbc, int ndims=NDIMS);
	dTensorBC5* clone(CopyMode::Enum copyMode)const
	{ return new dTensorBC5(*this, copyMode);}
        void copyfrom(const dTensorBC5& in);
	dTensorBC5(const dTensorBC5& in, CopyMode::Enum copyMode);
	int getmbc() const {return mbc;}
        int getsize(int n) const
        {
            switch(n)
            {
                case 5: return S5;
                case 4: return S4;
                case 3: return S3;
                case 2: return S2;
                case 1: return S1;
                default: return 1;
            }
        }
    private:
        int S1, S2, S3, S4, S5;
        int mbc;
        int ndims;
};

// Class with boundary conditions and four indices.
class dTensorBC4 : public dTensor4d
{

    public:
        dTensorBC4(int s1i, int s2i, int s3i, int s4i,
                int mbc, int ndims=NDIMS);
        dTensorBC4* clone(CopyMode::Enum copyMode)const
        { return new dTensorBC4(*this, copyMode);}
        void copyfrom(const dTensorBC4& in);
        int getmbc() const {return mbc;}
        int getsize(int n) const
        {
            switch(n)
            {
                case 4: return S4;
                case 3: return S3;
                case 2: return S2;
                case 1: return S1;
                default: return 1;
            }
        }

    private:

        // Methods:
        dTensorBC4(const dTensorBC4& in, CopyMode::Enum copyMode=CopyMode::DEEP);
        dTensorBC4& operator=(const dTensorBC4& in); // disabled

        // Fields:
        int S1, S2, S3, S4;
        int mbc;
        int ndims;

};

class dTensorBC3 : public dTensor3d
{
    private: // disabled
        dTensorBC3& operator=(const dTensorBC3& in);
    public:
        dTensorBC3(int s1i, int s2i, int s3i,
                int mbc, int ndims=NDIMS);
        void copyfrom(const dTensorBC3& in);
        int getmbc() const {return mbc;}
        int getsize(int n) const
        {
            switch(n)
            {
                case 3: return S3;
                case 2: return S2;
                case 1: return S1;
                default: return 1;
            }
        }
    private:
        int S1, S2, S3;
        int mbc;
        int ndims;
};

class dTensorBC2 : public dTensor2d
{
    private: // disabled
        dTensorBC2& operator=(const dTensorBC2& in);
    public:
        dTensorBC2(int s1i, int s2i,
                int mbc, int ndims=NDIMS);
        void copyfrom(const dTensorBC2& in);
        int getmbc() const {return mbc;}
        int getsize(int n) const
        {
            switch(n)
            {
                case 2: return S2;
                case 1: return S1;
                default: return 1;
            }
        }
    private:
        int S1, S2;
        int mbc;
        int ndims;
};

#endif
#ifndef tensors1d_h
#define tensors1d_h

#ifndef NDIMS
#include <dimdefs.h> // for NDIMS (number of spatial dimensions)
#endif
#if (!defined(NDIMS) || !(NDIMS>=0))
#error "NDIMS must be defined"
#endif

// What's the difference between DEEP and DIMS copy? (-DS)
namespace CopyMode
{
    enum Enum
    {
        DEEP = 1,
        DIMS = 2,
    };
}

class iTensorBase
{
    protected:
        void init();
        iTensorBase(){};
        iTensorBase(const iTensorBase& in);
    public:
        iTensorBase(int size_in):size(size_in){init();}
        ~iTensorBase();
        void setall(int);
        const int numel() const { return size; }
#ifdef CHECK_BOUNDS
        const int& vget(int k) const;
        int& vfetch(int k);
        void vset(int k, int value);
#else
        const int& vget(int k) const
        {
            return vec[k];
        }
        int& vfetch(int k) {
            return vec[k];
        }
        void vset(int k, int value){
            vec[k]=value;
        }
#endif
    protected:
        int* vec;
        int size;
};

class iTensor1d : public iTensorBase
{
    // data
    protected:
        int b1;
        // methods
    protected:
        iTensor1d(){};
    public:
        int getidx(int n1) const
        {
            int k = n1-b1;
            return k;
        }
        int getsize() const { return size; }
    public:
        // constructor takes size and initial index in each dimension
        iTensor1d( int s1i, int b1i ) : b1(b1i) { size=s1i; init(); }
        iTensor1d(const iTensor1d& in) : iTensorBase(in), b1(in.b1) {}

        const int& get(int n1) const
        { return vget(n1-b1); }
        void set(int n1, int value)
        { return vset(n1-b1, value); }
};

class iTensor1 : public iTensor1d
{
    public:
        iTensor1(int s1i) :
            iTensor1d(s1i,1) { }
};

// === major section: dtensors (arrays of double) ===

class dTensorBase
{
    private: // disabled
        dTensorBase& operator=(const dTensorBase& in);
    protected:
        void init();
        dTensorBase(){};
        dTensorBase(const dTensorBase& in, CopyMode::Enum copyMode=CopyMode::DEEP);
        void copyfrom(const dTensorBase& in);
    public:
        dTensorBase(int size_in):size(size_in){init();}
        ~dTensorBase();
        bool check();
        void setall(double);
        const int numel() const { return size; }
#ifdef CHECK_BOUNDS
        const double& vget(int k) const;
        double& vfetch(int k);
        void vset(int k, double value);
#else
        const double& vget(int k) const { return vec[k]; }
        double& vfetch(int k) {
            return vec[k];
#ifdef CHECK_INIT
            if(!(vec[k]==vec[k])) eprintf("vec[%d]=%24.16e",k,vec[k]);
#endif
        }
        void vset(int k, double value){vec[k]=value;}
#endif
    protected:
        double* vec;
        int size;
};

class dTensor1d : public dTensorBase
{
    // data
    protected:
        int b1;
        // methods
    private: // disabled
        dTensor1d& operator=(const dTensor1d& in);
    protected:
        dTensor1d(){};
    public:
        int getidx(int n1) const
        {
            int k = n1-b1;
            return k;
        }
        int getsize() const { return size; }
    public:
        // constructor takes size and initial index in each dimension
        dTensor1d( int s1i, int b1i ) : b1(b1i) { size=s1i; init(); }
        dTensor1d(const dTensor1d& in) : dTensorBase(in), b1(in.b1) {}
        void copyfrom(const dTensor1d& in);

        const double& get(int n1) const
        { return vget(n1-b1); }
        double& fetch(int n1)
        { return vfetch(n1-b1); }
        void set(int n1, double value)
        { return vset(n1-b1, value); }
};

class dTensor1 : public dTensor1d
{
    public:
        dTensor1(int s1i) :
            dTensor1d(s1i,1) { }
        void copyfrom(const dTensor1& in){ dTensor1d::copyfrom(in); }
    private: // disabled
        dTensor1& operator=(const dTensor1& in){ copyfrom(in); return *this; }
};

class dTensorBC1 : public dTensor1d
{
    private: // disabled
        dTensorBC1& operator=(const dTensorBC1& in);
    public:
        dTensorBC1(int s1i,
                int mbc, int ndims=1);
        void copyfrom(const dTensorBC1& in);
        int getmbc() const {return mbc;}
    private:
        int S1;
        int mbc;
        int ndims;
};

// methods for tensors
int count_fields(const char* str_in, char field_sep);
bool str_into_tensor(const char* str, iTensorBase& t, char field_sep);
bool str_into_tensor(const char* str, dTensorBase& t, char field_sep);
#endif
/// @file IniParams.cpp
/// Generated by /Users/seal/code/finess/python/finess/params/params.py

#include "IniParams.h"

#include "util.h"
#include "IniParser.h"

#include <fstream>
#include <string>

IniParams global_ini_params;


void IniParams::init(const std::string& inputFilename){
    using std::string;
    
    IniParser parser;
    {
        std::ifstream ifs(inputFilename.c_str());
	std::string ini_file_content = read_entire_stream(ifs);
	int parse_return_value = parser.parse(ini_file_content);
	if(parse_return_value != 0)
	    terminate("Error parsing " + inputFilename + ": line #" +
	              anyToString(parse_return_value));
    }

    this->ini_doc = parser.get_ini_doc();

// Defining code for member variables


    string global_alpha_str = this->ini_doc["finess"]["global_alpha"];
    
    if(global_alpha_str == ""){
        global_alpha_str = "false";
	this->ini_doc["finess"]["global_alpha"] = "false";
    }
        
    
    this->global_alpha = stringToAny<bool>(global_alpha_str);
        
    

    string mpp_limiter_str = this->ini_doc["finess"]["mpp_limiter"];
    
    if(mpp_limiter_str == ""){
        mpp_limiter_str = "false";
	this->ini_doc["finess"]["mpp_limiter"] = "false";
    }
        
    
    this->mpp_limiter = stringToAny<bool>(mpp_limiter_str);
        
    

    string output_dir_str = this->ini_doc["finess"]["output_dir"];
    
    if(output_dir_str == ""){
        output_dir_str = "output";
	this->ini_doc["finess"]["output_dir"] = "output";
    }
        
    
    this->output_dir = stringToAny<std::string>(output_dir_str);
        
    

    string ndims_str = this->ini_doc["finess"]["ndims"];
    
    if(ndims_str == "")
        terminate("finess.ndims is missing.");
    
    
    this->ndims = stringToAny<int>(ndims_str);
        
    

    string nout_str = this->ini_doc["finess"]["nout"];
    
    if(nout_str == ""){
        nout_str = "1";
	this->ini_doc["finess"]["nout"] = "1";
    }
        
    
    this->nout = stringToAny<int>(nout_str);
        
    

    string tfinal_str = this->ini_doc["finess"]["tfinal"];
    
    if(tfinal_str == "")
        terminate("finess.tfinal is missing.");
    
    
    this->tfinal = stringToAny<double>(tfinal_str);
        
    

    string initial_dt_str = this->ini_doc["finess"]["initial_dt"];
    
    if(initial_dt_str == "")
        terminate("finess.initial_dt is missing.");
    
    
    this->initial_dt = stringToAny<double>(initial_dt_str);
        
    

    string max_dt_str = this->ini_doc["finess"]["max_dt"];
    
    if(max_dt_str == "")
        terminate("finess.max_dt is missing.");
    
    
    this->max_dt = stringToAny<double>(max_dt_str);
        
    

    string desired_cfl_str = this->ini_doc["finess"]["desired_cfl"];
    
    if(desired_cfl_str == "")
        terminate("finess.desired_cfl is missing.");
    
    
    this->desired_cfl = stringToAny<double>(desired_cfl_str);
        
    

    string max_cfl_str = this->ini_doc["finess"]["max_cfl"];
    
    if(max_cfl_str == "")
        terminate("finess.max_cfl is missing.");
    
    
    this->max_cfl = stringToAny<double>(max_cfl_str);
        
    

    string nv_str = this->ini_doc["finess"]["nv"];
    
    if(nv_str == "")
        terminate("finess.nv is missing.");
    
    
    this->nv = stringToAny<int>(nv_str);
        
    

    string time_stepping_method_str = this->ini_doc["finess"]["time_stepping_method"];
    
    if(time_stepping_method_str == "")
        terminate("finess.time_stepping_method is missing.");
    
    
    this->time_stepping_method = stringToAny<IniParams::TimeSteppingMethod::enum_type>(time_stepping_method_str);
        
    
    if(this->time_stepping_method == IniParams::TimeSteppingMethod::DEFAULT)
        terminate("finess.time_stepping_method should be one of the following: SDC, Runge-Kutta, Lax-Wendroff, User-Defined, Multiderivative.");
               

    string space_order_str = this->ini_doc["finess"]["space_order"];
    
    if(space_order_str == "")
        terminate("finess.space_order is missing.");
    
    
    this->space_order = stringToAny<int>(space_order_str);
        
    

    string time_order_str = this->ini_doc["finess"]["time_order"];
    
    if(time_order_str == "")
        terminate("finess.time_order is missing.");
    
    
    this->time_order = stringToAny<int>(time_order_str);
        
    

    string verbosity_str = this->ini_doc["finess"]["verbosity"];
    
    if(verbosity_str == "")
        terminate("finess.verbosity is missing.");
    
    
    this->verbosity = stringToAny<int>(verbosity_str);
        
    

    string mcapa_str = this->ini_doc["finess"]["mcapa"];
    
    if(mcapa_str == "")
        terminate("finess.mcapa is missing.");
    
    
    this->mcapa = stringToAny<int>(mcapa_str);
        
    

    string maux_str = this->ini_doc["finess"]["maux"];
    
    if(maux_str == "")
        terminate("finess.maux is missing.");
    
    
    this->maux = stringToAny<int>(maux_str);
        
    

    string source_term_str = this->ini_doc["finess"]["source_term"];
    
    if(source_term_str == "")
        terminate("finess.source_term is missing.");
    
    
    this->source_term = stringToAny<bool>(source_term_str);
        
    

    string meqn_str = this->ini_doc["finess"]["meqn"];
    
    if(meqn_str == "")
        terminate("finess.meqn is missing.");
    
    
    this->meqn = stringToAny<int>(meqn_str);
        
    

    string weno_version_str = this->ini_doc["weno"]["weno_version"];
    
    if(weno_version_str == ""){
        weno_version_str = "JS";
	this->ini_doc["weno"]["weno_version"] = "JS";
    }
        
    
    this->weno_version = stringToAny<IniParams::WenoVersion::enum_type>(weno_version_str);
        
    
    if(this->weno_version == IniParams::WenoVersion::DEFAULT)
        terminate("weno.weno_version should be one of the following: Z, FD, JS.");
               

    string power_param_str = this->ini_doc["weno"]["power_param"];
    
    if(power_param_str == ""){
        power_param_str = "2.0";
	this->ini_doc["weno"]["power_param"] = "2.0";
    }
        
    
    this->power_param = stringToAny<double>(power_param_str);
        
    

    string alpha_scaling_str = this->ini_doc["weno"]["alpha_scaling"];
    
    if(alpha_scaling_str == ""){
        alpha_scaling_str = "1.1";
	this->ini_doc["weno"]["alpha_scaling"] = "1.1";
    }
        
    
    this->alpha_scaling = stringToAny<double>(alpha_scaling_str);
        
    

    string epsilon_str = this->ini_doc["weno"]["epsilon"];
    
    if(epsilon_str == ""){
        epsilon_str = "1e-06";
	this->ini_doc["weno"]["epsilon"] = "1e-06";
    }
        
    
    this->epsilon = stringToAny<double>(epsilon_str);
        
    

    string mx_str = this->ini_doc["grid"]["mx"];
    
    if(mx_str == "")
        terminate("grid.mx is missing.");
    
    
    this->mx = stringToAny<int>(mx_str);
        
    

    string my_str = this->ini_doc["grid"]["my"];
    
    if(my_str == "")
        terminate("grid.my is missing.");
    
    
    this->my = stringToAny<int>(my_str);
        
    

    string mbc_str = this->ini_doc["grid"]["mbc"];
    
    if(mbc_str == "")
        terminate("grid.mbc is missing.");
    
    
    this->mbc = stringToAny<int>(mbc_str);
        
    

    string xlow_str = this->ini_doc["grid"]["xlow"];
    
    if(xlow_str == "")
        terminate("grid.xlow is missing.");
    
    
    this->xlow = stringToAny<double>(xlow_str);
        
    

    string xhigh_str = this->ini_doc["grid"]["xhigh"];
    
    if(xhigh_str == "")
        terminate("grid.xhigh is missing.");
    
    
    this->xhigh = stringToAny<double>(xhigh_str);
        
    

    string ylow_str = this->ini_doc["grid"]["ylow"];
    
    if(ylow_str == "")
        terminate("grid.ylow is missing.");
    
    
    this->ylow = stringToAny<double>(ylow_str);
        
    

    string yhigh_str = this->ini_doc["grid"]["yhigh"];
    
    if(yhigh_str == "")
        terminate("grid.yhigh is missing.");
    
    
    this->yhigh = stringToAny<double>(yhigh_str);
        
    
    this->dx = (this->xhigh - this->xlow) / this->mx;
    this->dy = (this->yhigh - this->ylow) / this->my;

    string gamma_str = this->ini_doc["mhd"]["gamma"];
    
    if(gamma_str == "")
        terminate("mhd.gamma is missing.");
    
    
    this->gamma = stringToAny<double>(gamma_str);
        
    

    string rhol_str = this->ini_doc["initial"]["rhol"];
    
    if(rhol_str == "")
        terminate("initial.rhol is missing.");
    
    
    this->rhol = stringToAny<double>(rhol_str);
        
    

    string unl_str = this->ini_doc["initial"]["unl"];
    
    if(unl_str == "")
        terminate("initial.unl is missing.");
    
    
    this->unl = stringToAny<double>(unl_str);
        
    

    string utl_str = this->ini_doc["initial"]["utl"];
    
    if(utl_str == "")
        terminate("initial.utl is missing.");
    
    
    this->utl = stringToAny<double>(utl_str);
        
    

    string u3l_str = this->ini_doc["initial"]["u3l"];
    
    if(u3l_str == "")
        terminate("initial.u3l is missing.");
    
    
    this->u3l = stringToAny<double>(u3l_str);
        
    

    string pl_str = this->ini_doc["initial"]["pl"];
    
    if(pl_str == "")
        terminate("initial.pl is missing.");
    
    
    this->pl = stringToAny<double>(pl_str);
        
    

    string Bnl_str = this->ini_doc["initial"]["Bnl"];
    
    if(Bnl_str == "")
        terminate("initial.Bnl is missing.");
    
    
    this->Bnl = stringToAny<double>(Bnl_str);
        
    

    string Btl_str = this->ini_doc["initial"]["Btl"];
    
    if(Btl_str == "")
        terminate("initial.Btl is missing.");
    
    
    this->Btl = stringToAny<double>(Btl_str);
        
    

    string B3l_str = this->ini_doc["initial"]["B3l"];
    
    if(B3l_str == "")
        terminate("initial.B3l is missing.");
    
    
    this->B3l = stringToAny<double>(B3l_str);
        
    

    string rhor_str = this->ini_doc["initial"]["rhor"];
    
    if(rhor_str == "")
        terminate("initial.rhor is missing.");
    
    
    this->rhor = stringToAny<double>(rhor_str);
        
    

    string unr_str = this->ini_doc["initial"]["unr"];
    
    if(unr_str == "")
        terminate("initial.unr is missing.");
    
    
    this->unr = stringToAny<double>(unr_str);
        
    

    string utr_str = this->ini_doc["initial"]["utr"];
    
    if(utr_str == "")
        terminate("initial.utr is missing.");
    
    
    this->utr = stringToAny<double>(utr_str);
        
    

    string u3r_str = this->ini_doc["initial"]["u3r"];
    
    if(u3r_str == "")
        terminate("initial.u3r is missing.");
    
    
    this->u3r = stringToAny<double>(u3r_str);
        
    

    string pr_str = this->ini_doc["initial"]["pr"];
    
    if(pr_str == "")
        terminate("initial.pr is missing.");
    
    
    this->pr = stringToAny<double>(pr_str);
        
    

    string Bnr_str = this->ini_doc["initial"]["Bnr"];
    
    if(Bnr_str == "")
        terminate("initial.Bnr is missing.");
    
    
    this->Bnr = stringToAny<double>(Bnr_str);
        
    

    string Btr_str = this->ini_doc["initial"]["Btr"];
    
    if(Btr_str == "")
        terminate("initial.Btr is missing.");
    
    
    this->Btr = stringToAny<double>(Btr_str);
        
    

    string B3r_str = this->ini_doc["initial"]["B3r"];
    
    if(B3r_str == "")
        terminate("initial.B3r is missing.");
    
    
    this->B3r = stringToAny<double>(B3r_str);
        
    
// Checks

    if(ndims != 1 &&
       ndims != 2 &&
       ndims != 3)
            terminate("ndims is not one of [1, 2, 3].");


    if(!(this->nout > 0))
        terminate("nout should be > 0.");


    if(!(this->tfinal >= 0.0))
        terminate("tfinal should be >= 0.0.");


    if(!(this->initial_dt > 0.0))
        terminate("initial_dt should be > 0.0.");

if(!(this->max_dt >= this->initial_dt))
        terminate("finess.max_dt should > finess.initial_dt");
    

    if(!(this->desired_cfl > 0.0))
        terminate("desired_cfl should be > 0.0.");

if(!(this->max_cfl >= this->desired_cfl))
        terminate("finess.max_cfl should > this->desired_cfl");
    

    if(!(this->nv > 0))
        terminate("nv should be > 0.");

    if(space_order != 1 &&
       space_order != 3 &&
       space_order != 5 &&
       space_order != 7 &&
       space_order != 9 &&
       space_order != 11)
            terminate("space_order is not one of [1, 3, 5, 7, 9, 11].");

    if(time_order != 1 &&
       time_order != 2 &&
       time_order != 3 &&
       time_order != 4 &&
       time_order != 5)
            terminate("time_order is not one of [1, 2, 3, 4, 5].");

    if(verbosity != 0 &&
       verbosity != 1)
            terminate("verbosity is not one of [0, 1].");


    if(!(this->mcapa >= 0))
        terminate("mcapa should be >= 0.");

if(!(this->maux >= this->mcapa))
        terminate("finess.maux should >= finess.mcapa");
    

    if(!(this->meqn >= 1))
        terminate("meqn should be >= 1.");


    if(!(this->alpha_scaling >= 1.0))
        terminate("alpha_scaling should be >= 1.0.");


    if(!(this->epsilon > 0.0))
        terminate("epsilon should be > 0.0.");


    if(!(this->mx > 0))
        terminate("mx should be > 0.");


    if(!(this->my > 0))
        terminate("my should be > 0.");


    if(!(this->mbc >= 0))
        terminate("mbc should be >= 0.");

if(!(xhigh > xlow))
        terminate("grid.xhigh should > grid.xlow.");
if(!(yhigh > ylow))
        terminate("grid.yhigh should > grid.ylow.");

    if(!(this->gamma > 0.0))
        terminate("gamma should be > 0.0.");


    if(!(this->rhol > 0.0))
        terminate("rhol should be > 0.0.");


    if(!(this->pl > 0.0))
        terminate("pl should be > 0.0.");


    if(!(this->rhor > 0.0))
        terminate("rhor should be > 0.0.");


    if(!(this->pr > 0.0))
        terminate("pr should be > 0.0.");


}

#ifndef _INIPARAMS_H_
#define _INIPARAMS_H_
/// @file IniParams.h
/// Generated by /Users/seal/code/finess/python/finess/params/params.py
    

#include <string>
#include "util.h"
#include "IniParser.h"

class IniParams;
extern IniParams global_ini_params;


class IniParams{
public:
    void init(const std::string& inputFilename);
private:
    IniParser::ini_doc_type ini_doc;
public:
    std::string ini_doc_as_string(){
        return IniParser::ini_doc_as_string(this->ini_doc);
    }
// Type definitions 
public:
        struct WenoVersion{
            enum enum_type {Z, FD, JS, DEFAULT};
        };
        
public:
        struct TimeSteppingMethod{
            enum enum_type {SDC, RK, LxW, USER_DEFINED, MD, DEFAULT};
        };
        

// Member variables declarations
private:
        bool global_alpha;
        
private:
        bool mpp_limiter;
        
private:
        std::string output_dir;
        
private:
        int ndims;
        
private:
        int nout;
        
private:
        double tfinal;
        
private:
        double initial_dt;
        
private:
        double max_dt;
        
private:
        double desired_cfl;
        
private:
        double max_cfl;
        
private:
        int nv;
        
private:
        IniParams::TimeSteppingMethod::enum_type time_stepping_method;
        
private:
        int space_order;
        
private:
        int time_order;
        
private:
        int verbosity;
        
private:
        int mcapa;
        
private:
        int maux;
        
private:
        bool source_term;
        
private:
        int meqn;
        
private:
        IniParams::WenoVersion::enum_type weno_version;
        
private:
        double power_param;
        
private:
        double alpha_scaling;
        
private:
        double epsilon;
        
private:
        int mx;
        
private:
        int my;
        
private:
        int mbc;
        
private:
        double xlow;
        
private:
        double xhigh;
        
private:
        double ylow;
        
private:
        double yhigh;
        
private:
        double dx;
        
private:
        double dy;
        
private:
        double gamma;
        
private:
        double rhol;
        
private:
        double unl;
        
private:
        double utl;
        
private:
        double u3l;
        
private:
        double pl;
        
private:
        double Bnl;
        
private:
        double Btl;
        
private:
        double B3l;
        
private:
        double rhor;
        
private:
        double unr;
        
private:
        double utr;
        
private:
        double u3r;
        
private:
        double pr;
        
private:
        double Bnr;
        
private:
        double Btr;
        
private:
        double B3r;
        

// Accessor definitions
public:
        inline bool get_global_alpha(){
            return this->global_alpha;
        }
        
public:
        inline bool get_mpp_limiter(){
            return this->mpp_limiter;
        }
        
public:
        inline std::string get_output_dir(){
            return this->output_dir;
        }
        
public:
        inline int get_ndims(){
            return this->ndims;
        }
        
public:
        inline int get_nout(){
            return this->nout;
        }
        
public:
        inline double get_tfinal(){
            return this->tfinal;
        }
        
public:
        inline double get_initial_dt(){
            return this->initial_dt;
        }
        
public:
        inline double get_max_dt(){
            return this->max_dt;
        }
        
public:
        inline double get_desired_cfl(){
            return this->desired_cfl;
        }
        
public:
        inline double get_max_cfl(){
            return this->max_cfl;
        }
        
public:
        inline int get_nv(){
            return this->nv;
        }
        
public:
        inline IniParams::TimeSteppingMethod::enum_type get_time_stepping_method(){
            return this->time_stepping_method;
        }
        
public:
        inline int get_space_order(){
            return this->space_order;
        }
        
public:
        inline int get_time_order(){
            return this->time_order;
        }
        
public:
        inline int get_verbosity(){
            return this->verbosity;
        }
        
public:
        inline int get_mcapa(){
            return this->mcapa;
        }
        
public:
        inline int get_maux(){
            return this->maux;
        }
        
public:
        inline bool get_source_term(){
            return this->source_term;
        }
        
public:
        inline int get_meqn(){
            return this->meqn;
        }
        
public:
        inline IniParams::WenoVersion::enum_type get_weno_version(){
            return this->weno_version;
        }
        
public:
        inline double get_power_param(){
            return this->power_param;
        }
        
public:
        inline double get_alpha_scaling(){
            return this->alpha_scaling;
        }
        
public:
        inline double get_epsilon(){
            return this->epsilon;
        }
        
public:
        inline int get_mx(){
            return this->mx;
        }
        
public:
        inline int get_my(){
            return this->my;
        }
        
public:
        inline int get_mbc(){
            return this->mbc;
        }
        
public:
        inline double get_xlow(){
            return this->xlow;
        }
        
public:
        inline double get_xhigh(){
            return this->xhigh;
        }
        
public:
        inline double get_ylow(){
            return this->ylow;
        }
        
public:
        inline double get_yhigh(){
            return this->yhigh;
        }
        
public:
        inline double get_dx(){
            return this->dx;
        }
        
public:
        inline double get_dy(){
            return this->dy;
        }
        
public:
        inline double get_gamma(){
            return this->gamma;
        }
        
public:
        inline double get_rhol(){
            return this->rhol;
        }
        
public:
        inline double get_unl(){
            return this->unl;
        }
        
public:
        inline double get_utl(){
            return this->utl;
        }
        
public:
        inline double get_u3l(){
            return this->u3l;
        }
        
public:
        inline double get_pl(){
            return this->pl;
        }
        
public:
        inline double get_Bnl(){
            return this->Bnl;
        }
        
public:
        inline double get_Btl(){
            return this->Btl;
        }
        
public:
        inline double get_B3l(){
            return this->B3l;
        }
        
public:
        inline double get_rhor(){
            return this->rhor;
        }
        
public:
        inline double get_unr(){
            return this->unr;
        }
        
public:
        inline double get_utr(){
            return this->utr;
        }
        
public:
        inline double get_u3r(){
            return this->u3r;
        }
        
public:
        inline double get_pr(){
            return this->pr;
        }
        
public:
        inline double get_Bnr(){
            return this->Bnr;
        }
        
public:
        inline double get_Btr(){
            return this->Btr;
        }
        
public:
        inline double get_B3r(){
            return this->B3r;
        }
        
};
//stringToAny specializations
template<>
inline IniParams::WenoVersion::enum_type stringToAny<IniParams::WenoVersion::enum_type>(const std::string& s){
    return
        s == "Z" ? IniParams::WenoVersion::Z :
        s == "FD" ? IniParams::WenoVersion::FD :
        s == "JS" ? IniParams::WenoVersion::JS :
        IniParams::WenoVersion::DEFAULT;
        
}
        
template<>
inline IniParams::TimeSteppingMethod::enum_type stringToAny<IniParams::TimeSteppingMethod::enum_type>(const std::string& s){
    return
        s == "SDC" ? IniParams::TimeSteppingMethod::SDC :
        s == "Runge-Kutta" ? IniParams::TimeSteppingMethod::RK :
        s == "Lax-Wendroff" ? IniParams::TimeSteppingMethod::LxW :
        s == "User-Defined" ? IniParams::TimeSteppingMethod::USER_DEFINED :
        s == "Multiderivative" ? IniParams::TimeSteppingMethod::MD :
        IniParams::TimeSteppingMethod::DEFAULT;
        
}
        
#endif
#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

    const double gamma1 = global_ini_params.get_gamma() - 1.0;


    const double rhol = global_ini_params.get_rhol();
    const double u1l  = global_ini_params.get_unl();
    const double u2l  = global_ini_params.get_utl();
    const double u3l  = global_ini_params.get_u3l();
    const double pl   = global_ini_params.get_pl();
    const double B1l  = global_ini_params.get_Bnl();
    const double B2l  = global_ini_params.get_Btl();
    const double B3l  = global_ini_params.get_B3l();

    const double rhor = global_ini_params.get_rhor();
    const double u1r  = global_ini_params.get_unr();
    const double u2r  = global_ini_params.get_utr();
    const double u3r  = global_ini_params.get_u3r();
    const double pr   = global_ini_params.get_pr();
    const double B1r  = global_ini_params.get_Bnr();
    const double B2r  = global_ini_params.get_Btr();
    const double B3r  = global_ini_params.get_B3r();

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i,1);
      double rho,u1,u2,u3,press,B1,B2,B3,energy;

      if(x<0.0e0)
        {
	  rho   =  rhol;
	  u1    =  u1l;
	  u2    =  u2l;
	  u3    =  u3l;
	  press =  pl;
	  B1    =  B1l;
	  B2    =  B2l;
	  B3    =  B3l;
        }
      else
        {
	  rho   =  rhor;
	  u1    =  u1r;
	  u2    =  u2r;
	  u3    =  u3r;
	  press =  pr;
	  B1    =  B1r;
	  B2    =  B2r;
	  B3    =  B3r;
        }

      energy = press/gamma1 
	+ 0.5*rho*(u1*u1 + u2*u2 + u3*u3)
	+ 0.5*(B1*B1 + B2*B2 + B3*B3);

      qvals.set(i,1, rho );     // density
      qvals.set(i,2, rho*u1 );  // 1-momentum
      qvals.set(i,3, rho*u2 );  // 2-momentum
      qvals.set(i,4, rho*u3 );  // 3-momentum
      qvals.set(i,5, energy );  // energy
      qvals.set(i,6, B1 );      // B1
      qvals.set(i,7, B2 );      // B2
      qvals.set(i,8, B3 );      // B3      
    }
}
#include "tensors.h"
#include "StateVars.h"

// This is a user-supplied routine that sets the the boundary conditions
//
// TODO - what type of boundary conditions are being applied here? -DS
//
// See also: ...
void SetBndValues( StateVars& Q )
{

    dTensorBC3& q     = Q.ref_q();
    dTensorBC3& aux   = Q.ref_aux();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for(int i=0; i>=(1-mbc); i--)
    {
        for(int j=1; j<=my; j++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(1,j,m);
                q.set(i,j,m, tmp );
            }

            // aux values
            for(int m=1; m<=maux; m++)
            {
                double tmp = aux.get(1,j,m);
                aux.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************



    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for(int i=(mx+1); i<=(mx+mbc); i++)
    {
        for(int j=1; j<=my; j++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(mx,j,m);
                q.set(i,j,m, tmp );
            }

            // aux values
            for(int m=1; m<=maux; m++)
            {
                double tmp = aux.get(mx,j,m);
                aux.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************



    // ***********************************************
    // BOTTOM BOUNDARY
    // ***********************************************
    for(int j=0; j>=(1-mbc); j--)
    {
        for(int i=1; i<=mx; i++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(i,1,m);
                q.set(i,j,m, tmp );
            }               

            // aux values
            for(int m=1; m<=maux; m++)
            {
                double tmp = aux.get(i,1,m);
                aux.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************



    // ***********************************************
    // TOP BOUNDARY
    // ***********************************************
    for(int j=(my+1); j<=(my+mbc); j++)
    {
        for(int i=1; i<=mx; i++)
        {           
            // q values
            for(int m=1; m<=meqn; m++)
            {
                double tmp = q.get(i,my,m);
                q.set(i,j,m, tmp );
            }

            // aux values
            for(int m=1; m<=maux; m++)
            {
                double tmp = aux.get(i,my,m);
                aux.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************


    // ***********************************************
    // BOTTOM LEFT CORNER
    // ***********************************************
    for(int i=1; i<=mbc; i++)
        for(int j=1; j<=mbc; j++)
        {
            for(int m=1; m<=meqn; m++)
            {     
                q.set(1-i,1-j,m, q.get(1,1,m) );
            }
            for(int m=1; m<=maux; m++)
            {     
                aux.set(1-i,1-j,m, aux.get(1,1,m) );
            }
        }
    // ***********************************************


    // ***********************************************
    // BOTTOM RIGHT CORNER
    // ***********************************************
    for(int i=1; i<=mbc; i++)
        for(int j=1; j<=mbc; j++)
        {
            for(int m=1; m<=meqn; m++)
            {     
                q.set(mx+i,1-j,m, q.get(mx,1,m) );
            }
            for(int m=1; m<=maux; m++)
            {     
                aux.set(mx+i,1-j,m, aux.get(mx,1,m) );
            }
        }
    // ***********************************************


    // ***********************************************
    // TOP RIGHT CORNER
    // ***********************************************
    for(int i=1; i<=mbc; i++)
        for(int j=1; j<=mbc; j++)
        {
            for(int m=1; m<=meqn; m++)
            {     
                q.set(mx+i,my+j,m, q.get(mx,my,m) );
            }
            for(int m=1; m<=maux; m++)
            {     
                aux.set(mx+i,my+j,m, aux.get(mx,my,m) );
            }
        }
    // ***********************************************


    // ***********************************************
    // TOP LEFT CORNER
    // ***********************************************
    for(int i=1; i<=mbc; i++)
        for(int j=1; j<=mbc; j++)
        {
            for(int m=1; m<=meqn; m++)
            {     
                q.set(1-i,my+j,m, q.get(1,my,m) );
            }
            for(int m=1; m<=maux; m++)
            {     
                aux.set(1-i,my+j,m, aux.get(1,my,m) );
            }
        }
    // ***********************************************

}

void SetBndValuesX( StateVars& Q )
{ SetBndValues( Q ); }

void SetBndValuesY( StateVars& Q )
{ SetBndValues( Q ); }
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals,
		    const dTensor2& auxvals, dTensor2& source)
{
  int i,m;
  int numpts=xpts.getsize(1);
  int meqn=qvals.getsize(2);
  double x,y;
  
  for (i=1; i<=numpts; i++)
    {
      x = xpts.get(i,1);
      y = xpts.get(i,2);
        
      for (m=1; m<=meqn; m++)
	{
	  source.set(i,m, 0.0e0 );
	}
    }
}
// =========================================================================
//
//  --------------------------------------------
//  FINESS: See LICENSE.txt for licensing details
//  --------------------------------------------
//
//    Lead Developer:  
//             David C. Seal
//             Michigan University
//             Department of Mathematics
//             619 Red Cedar Road
//             East Lansing, MI 48823
//             seal@math.msu.edu
//
// =========================================================================

int main(int argc, char* argv[])
{

    // NOTE: You should not have to modify this part of the code.
    //
    //       To change parameters, modify the following files:
    //            1. parameters.ini -- basic data file, can modify
    //                  number of grid points, time step, order of
    //                  accuracy in both space and time, etc...
    //            2. QinitFunc.cpp          -- initial condition file
    //            3. AuxFunc.cpp            -- auxiliary variable file
    //            4. SourceTermFunc.pp      -- source term file
    //            5. FluxFunc.cpp           -- flux function file
    //            6. SetWaveSpd.cpp         -- eigenvalues of flux Jacobian file
    //            7. ProjectLeftEig.cpp     -- left eigenvectors of flux Jacobian file
    //            8. ProjectRightEig.cpp    -- right eigenvectors of flux Jacobian file
    //            9. SetBndValues.cpp       -- boundary conditions files
    //

    int main_global(int argc, char* argv[]);
    return main_global(argc,argv);

}
