#include <cmath>
#include <stdexcept>
#include "assert.h"
#include "tensors.h"
#include "CentralDifferences.h"
#include "IniParams.h"

central_differences_t GetCentralDifferences()
{

    // TODO - implement WENO3 and WENO-Z version
    if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 5)
        return &WenoReconstruct_JS5;
//  else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::JS && global_ini_params.get_space_order() == 7)
//      return &WenoReconstruct_JS7;
//  else if(global_ini_params.get_weno_version() == IniParams::WenoVersion::Z  && global_ini_params.get_space_order() == 9)
//  {
//      printf("Warning: we're not sure these are the correct coefficients for WENOZ-9!\n");
//      return &WenoReconstruct_Z9;
//  }
    else
        throw(std::logic_error("Requested Central Difference Operator is not implemented."));
}

// ------------------------------------------------- // 
// SECTION: Central Finite difference approximations //
// ------------------------------------------------- // 

// Compute all derivatives that come from a 5 point stencil
//
// Input:
//
//      dx                 - cell width
//
//      f( 1:mcomps, 1:ws ) - list of meqn functions to be differentiated. ws =
//                          size of stencil under consideration.  ws = 5 for
//                          now.
//
// Output:
//
//      fderivs( 1:mcomps, 1:ws ) - The differentiated function.
//
//                                  fderivs(:,1) = f,
//                                  fderivs(:,2) = fx,
//                                  fderivs(:,3) = fxx,
//                                  fderivs(:,4) = fxxx,
//                                  fderivs(:,5) = fxxxx.
//
void CentralDifferences5( double dx, const dTensor2& f, dTensor2& fderivs )
{

    const int mcomps = f.getsize( 1 );  // Usually "meqn"
    const int STENCIL_SIZE = 5;
    assert_eq( f.getsize( 2 ), STENCIL_SIZE );

    for( int m=1; m <= mcomps; m++ )
    {

        // Zeroth derivative
        fderivs.set(m,1, f.get(m, (STENCIL_SIZE/2)+1 ) );

        // First and higher-derivatives
        for( int nderiv=1; nderiv <= STENCIL_SIZE-1; nderiv++ )
        {
            double tmp = 0.;
            for( int i=1; i <= STENCIL_SIZE; i++ )
            {
                tmp += f.get( m, i ) * deriv_matrix5[nderiv-1][i-1];
            }

            fderivs.set( m, nderiv+1, tmp*pow(dx,-nderiv) );
        }
    }

}

// Differences with a 7-point stencil
void CentralDifferences7( double dx, const dTensor2& f, dTensor2& fderivs )
{

    const int mcomps = f.getsize( 1 );  // Usually "meqn"
    const int STENCIL_SIZE = 7;

    assert_eq( f.getsize( 2 ), STENCIL_SIZE );
    for( int m=1; m <= mcomps; m++ )
    {

        // Zeroth derivative
        fderivs.set(m,1, f.get(m, (STENCIL_SIZE/2)+1 ) );

        // First and higher-derivatives
        for( int nderiv=1; nderiv <= STENCIL_SIZE-1; nderiv++ )
        {
            double tmp = 0.;
            for( int i=1; i <= STENCIL_SIZE; i++ )
            {
                tmp += f.get( m, i ) * deriv_matrix7[nderiv-1][i-1];
            }

            fderivs.set( m, nderiv+1, tmp*pow(dx,-nderiv) );
        }
    }

}

// Differences with a 9-point stencil
void CentralDifferences9( double dx, const dTensor2& f, dTensor2& fderivs )
{

    const int mcomps = f.getsize( 1 );  // Usually "meqn" but doesn't have to be
    const int STENCIL_SIZE = 9;

    assert_eq( f.getsize( 2 ), STENCIL_SIZE );
    for( int m=1; m <= mcomps; m++ )
    {

        // Zeroth derivative
        fderivs.set(m,1, f.get(m, (STENCIL_SIZE/2)+1 ) );

        // First and higher-derivatives
        for( int nderiv=1; nderiv <= STENCIL_SIZE-1; nderiv++ )
        {
            double tmp = 0.;
            for( int i=1; i <= STENCIL_SIZE; i++ )
            {
                tmp += f.get( m, i ) * deriv_matrix9[nderiv-1][i-1];
            }

            fderivs.set( m, nderiv+1, tmp*pow(dx,-nderiv) );
        }
    }

}
