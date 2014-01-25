#include <iostream>
#include "dog_math.h"
#include "tensors.h"
#include "stdlib.h"
#include "DogParamsCart2.h"

using namespace std;

double GetCFL(double dt, double dtmax,
        const dTensor2& prim_vol,
        const dTensorBC3& aux,
        const dTensorBC1& smax)
{

    const int mx = dogParamsCart2.get_mx();
    const int my = dogParamsCart2.get_my();
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
        cfl = Max( dt*smax.get(j) / prim_vol.get(i,j), cfl);  
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
