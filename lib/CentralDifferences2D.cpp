#include <cmath>
#include <stdexcept>
#include "assert.h"
#include "tensors.h"
#include "CentralDifferences2D.h"
#include "IniParams.h"
#include "dog_math.h"

central_differences_t2d GetCentralDifferences2D()
{

    // TODO - implement WENO3 and WENO-Z version
    if( global_ini_params.get_space_order() == 5)
        return &CentralDifferences2D5;
    else if( global_ini_params.get_space_order() == 7)
        return &CentralDifferences2D7;
    else if( global_ini_params.get_space_order() == 9)
        return &CentralDifferences2D9;
    else if( global_ini_params.get_space_order() == 11)
        return &CentralDifferences2D11;
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
/*
[0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4]
[0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4]
[f f_x f_xx f_3x f_4x f_y f_xy f_2xy f_3xy f_4xy f_2y f_x2y f_2x2y f_3x2y f_4x2y f_3y f_x3y f_2x3y f_3x3y f_4x3y f_4y f_x4y f_2x4y f_3x4y f_4x4y]
*/

//
void CentralDifferences2D5( double dx,double dy, const dTensor3& f, dTensor3& fderivs )
{

    const int mcomps = f.getsize( 1 );  // Usually "meqn"
    const int STENCIL_SIZE = 5;
    assert_eq( f.getsize( 2 ), STENCIL_SIZE );
    assert_eq( f.getsize( 3 ), STENCIL_SIZE );

    int nx=(1.0/dx-1.0);
    int ny=(1.0/dy-1.0);


    for( int m=1; m <= mcomps; m++ )
    {

        // First and higher-derivatives
        for( int nderivx=1; nderivx <= STENCIL_SIZE; nderivx++ )
        for( int nderivy=1; nderivy <= STENCIL_SIZE; nderivy++ )
        {
            double tmp = 0.;
            int nderiv = ((nderivy-1)*STENCIL_SIZE+nderivx);
            for( int j=1; j <= STENCIL_SIZE; j++ )
            {
              for(int i=1; i<=STENCIL_SIZE;  i++)
               { tmp += f.get( m, i, j ) * Mdderiv_matrix5[nderiv-1][(j-1)*STENCIL_SIZE+(i-1)];}
            }
            //No longer need dx or dy as series is in xi eta
            fderivs.set( m, nderivx,nderivy,tmp/factorial[(int)Ix5[nderiv-1]]/factorial[(int)Iy5[nderiv-1]]);
        }
        
        // Zeroth derivative
        fderivs.set(m,1,1, f.get(m, (STENCIL_SIZE/2)+1, (STENCIL_SIZE/2)+1 ) );
    }

}

// Differences with a 7-point stencil
void CentralDifferences2D7( double dx,double dy, const dTensor3& f, dTensor3& fderivs )
{

    const int mcomps = f.getsize( 1 );  // Usually "meqn"
    const int STENCIL_SIZE = 7;
    int nx=(1.0/dx-1.0);
    int ny=(1.0/dy-1.0);


    assert_eq( f.getsize( 2 ), STENCIL_SIZE );
    assert_eq( f.getsize( 3 ), STENCIL_SIZE );

    for( int m=1; m <= mcomps; m++ )
    {

        // First and higher-derivatives
        for( int nderivx=1; nderivx <= STENCIL_SIZE; nderivx++ )
        for( int nderivy=1; nderivy <= STENCIL_SIZE; nderivy++ )
        {
            double tmp = 0.;
            int nderiv = ((nderivy-1)*STENCIL_SIZE+nderivx);
            for( int i=1; i <= STENCIL_SIZE; i++ )
            {
              for(int j=1; j<=STENCIL_SIZE;  j++)
               { tmp += (f.get( m, i, j ) ) * Mdderiv_matrix7[nderiv-1][(j-1)*STENCIL_SIZE+(i-1)];}
            }
            //No longer need dx or dy as series is in xi eta
            fderivs.set( m, nderivx,nderivy,tmp/factorial[(int)Ix7[nderiv-1]]/factorial[(int)Iy7[nderiv-1]]);
        }
        // Zeroth derivative
        fderivs.set(m,1,1, f.get(m, (STENCIL_SIZE/2)+1, (STENCIL_SIZE/2)+1 ) );
    }


}

// Differences with a 9-point stencil
void CentralDifferences2D9( double dx,double dy, const dTensor3& f, dTensor3& fderivs )
{

    const int mcomps = f.getsize( 1 );  // Usually "meqn" but doesn't have to be
    const int STENCIL_SIZE = 9;
    int nx=(1.0/dx-1.0);
    int ny=(1.0/dy-1.0);


    assert_eq( f.getsize( 2 ), STENCIL_SIZE );
    assert_eq( f.getsize( 3 ), STENCIL_SIZE );

    for( int m=1; m <= mcomps; m++ )
    {

        // First and higher-derivatives
        for( int nderivx=1; nderivx <= STENCIL_SIZE; nderivx++ )
        for( int nderivy=1; nderivy <= STENCIL_SIZE; nderivy++ )
        {
            double tmp = 0.;
            int nderiv = ((nderivy-1)*STENCIL_SIZE+nderivx);
            for( int i=1; i <= STENCIL_SIZE; i++ )
            {
              for(int j=1; j<=STENCIL_SIZE;  j++)
               { tmp += f.get( m, i, j ) * Mdderiv_matrix9[nderiv-1][(j-1)*STENCIL_SIZE+(i-1)];}
            }
            //No longer need dx or dy as series is in xi eta
            fderivs.set( m, nderivx,nderivy,tmp/factorial[(int)Ix9[nderiv-1]]/factorial[(int)Iy9[nderiv-1]]);
        }
        // Zeroth derivative
        fderivs.set(m,1,1, f.get(m, (STENCIL_SIZE/2)+1, (STENCIL_SIZE/2)+1 ) );
    }


}

// Differences with an 11-point stencil
void CentralDifferences2D11( double dx, double dy, const dTensor3& f, dTensor3& fderivs )
{

    const int mcomps = f.getsize( 1 );  // Usually "meqn" but doesn't have to be
    const int STENCIL_SIZE = 11;
    int nx=(1.0/dx-1.0);
    int ny=(1.0/dy-1.0);

    assert_eq( f.getsize( 2 ), STENCIL_SIZE );
    assert_eq( f.getsize( 3 ), STENCIL_SIZE );

    for( int m=1; m <= mcomps; m++ )
    {

        // First and higher-derivatives
        for( int nderivx=1; nderivx <= STENCIL_SIZE; nderivx++ )
        for( int nderivy=1; nderivy <= STENCIL_SIZE; nderivy++ )
        {
            double tmp = 0.;
            int nderiv = ((nderivy-1)*STENCIL_SIZE+nderivx);
            for( int i=1; i <= STENCIL_SIZE; i++ )
            {
              for(int j=1; j<=STENCIL_SIZE;  j++)
               { tmp += f.get( m, i, j ) * Mdderiv_matrix11[nderiv-1][(j-1)*STENCIL_SIZE+(i-1)];}
            }
            //No longer need dx or dy as series is in xi eta
            fderivs.set( m, nderivx,nderivy,tmp/factorial[(int)Ix11[nderiv-1]]/factorial[(int)Iy11[nderiv-1]]);

        }
        // Zeroth derivative
        fderivs.set(m,1,1, f.get(m, (STENCIL_SIZE/2)+1, (STENCIL_SIZE/2)+1 ) );
    }



}
