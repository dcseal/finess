#include <iostream>
#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"
#include "assert.h"
#include "StateVars.h"

using namespace std;

void FluxFunc(const dTensor2&,const dTensor2&,const dTensor2&,dTensor3&);

// Used for construcing the flux function
void SampleFunctionTypeB( 
        int istart, int iend,
        int jstart, int jend,
        int kstart, int kend,
        const dTensorBC4& qin, 
        const dTensorBC4& auxin,  
        dTensorBC5& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&));



// Same as above, but this one only returns flux values, and not just the RHS.
void ConstructLFL( const double dt, const StateVars& Qnew, 
    dTensorBC4& fLF, dTensorBC4& gLF, dTensorBC4& hLF, dTensorBC4& Lstar, dTensorBC4& smax )
{

    const dTensorBC4&    q = Qnew.const_ref_q  ();
    const dTensorBC4&  aux = Qnew.const_ref_aux();

//    void SetBndValues(StateVars& Qnew); 
//    SetBndValues( Qnew );

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int mz     = global_ini_params.get_mz();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double     dz = global_ini_params.get_dz();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();
    const double   zlow = global_ini_params.get_zlow();

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC5 R( mx, my, mz, meqn, 3, mbc );  // place-holder for the flux function

    const int mbc_small      = 3;
    void GlobalWaveSpd(
	    const dTensorBC4& q, 
	    const dTensorBC4& aux, 
	    double& alpha1, double& alpha2, double& alpha3 );

    // Compute finite difference approximations on all of the conserved
    // variables:
    SampleFunctionTypeB( 1-mbc, mx+mbc, 1-mbc, my+mbc, 1-mbc, mz+mbc, q, aux, R, &FluxFunc );

    // Find a global wave speed
    double alpha_x, alpha_y, alpha_z;
    GlobalWaveSpd(q, aux, alpha_x, alpha_y, alpha_z );
    smax.set(1,1,1, 1, alpha_x);
    smax.set(1,1,1, 2, alpha_y);
    smax.set(1,1,1, 3, alpha_z);

    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
    for (int k = 1; k <= mz;   k++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =  0.5*(R.get(i-1,j,k,m,1)+R.get(i,j,k,m,1)) 
                   + 0.5*alpha_x*(q.get(i-1,j,k,m)-q.get(i,j,k,m));
        fLF.set(i,j,k,m, hf );
    }

    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
    for (int k = 1; k<= mz;   k++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =  0.5*(R.get(i,j-1,k,m,2)+R.get(i,j,k,m,2)) 
                   + 0.5*alpha_y*(q.get(i,j-1,k,m)-q.get(i,j,k,m));
        gLF.set(i,j,k,m, hf );
    }

    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my;   j++)
    for (int k = 1; k<= mz+1; k++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =  0.5*(R.get(i,j,k-1,m,3)+R.get(i,j,k,m,3)) 
                   + 0.5*alpha_z*(q.get(i,j,k-1,m)-q.get(i,j,k,m));
        hLF.set(i,j,k,m, hf );
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
	for (int k=1; k<=mz; k++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(fLF.get(i+1,j,k, m) - fLF.get(i,j,k,m) ) / dx;
                tmp =  tmp   -(gLF.get(i,j+1,k, m) - gLF.get(i,j,k,m) ) / dy;
                tmp =  tmp   -(hLF.get(i,j,k+1, m) - hLF.get(i,j,k,m) ) / dz;
                Lstar.set(i,j,k,m, tmp );
            }
        }
    }

}
