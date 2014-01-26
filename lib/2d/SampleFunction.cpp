#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "stdlib.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
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

    const int meqn    = dogParams.get_meqn();
    const int maux    = dogParams.get_maux();
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
    const double dx   = dogParamsCart2.get_dx();
    const double dy   = dogParamsCart2.get_dy();
    const double xlow = dogParamsCart2.get_xlow();
    const double ylow = dogParamsCart2.get_ylow();

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

            // Evaluate integrals
            for (int m1=1; m1<=mlength; m1++)
            {
                Fout.set(i, j, m1, fvals.get(1,  m1) );
            }

        }
    }

}
