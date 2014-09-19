#include <iostream>
#include "dog_math.h"
#include "tensors.h"
#include "stdlib.h"
#include "IniParams.h"

using namespace std;

double GetCFL(double dt, double dtmax,
        const dTensorBC2& aux,
        const dTensorBC1& smax)
{

    const int    mx = global_ini_params.get_mx();
    const double dx = global_ini_params.get_dx();

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
    for (int j=1; j<=mx; j++)
    {
        cfl = Max( dt*smax.get(j) / dx, cfl);  
    }

    if (cfl>1.0e8)
    {
        cout << endl;
        cout << " Error: CFL number is out of bounds ... " << endl;
        cout << "     CFL = " << cfl << endl;
        cout << endl;
        exit(1);
    }

    return cfl;
}
