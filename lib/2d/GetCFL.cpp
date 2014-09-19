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
