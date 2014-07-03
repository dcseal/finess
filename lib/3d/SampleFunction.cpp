#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "stdlib.h"
#include "DogParams.h"
#include "DogParamsCart3.h"
using namespace std;

// All-purpose routine for sampling a function, and saving its data into a
// single tensor.
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//
// TODO - change the indexing here!
//
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
    int kstart, int kend,
    const dTensorBC4& qin, 
    const dTensorBC4& auxin,  
          dTensorBC4& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&))
{    

    const int meqn    = dogParams.get_meqn();
    const int maux    = dogParams.get_maux();
    const int mlength = Fout.getsize(4);

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
        cout << "       jstart = " << jstart << endl;
        cout << "         jend = " << jend << endl;
        cout << "       kstart = " << kstart << endl;
        cout << "         kend = " << kend << endl;
        cout << endl;
        exit(1);
    }

    // Grid information
    const double dx   = dogParamsCart3.get_dx();
    const double dy   = dogParamsCart3.get_dy();
    const double dz   = dogParamsCart3.get_dz();
    const double xlow = dogParamsCart3.get_xlow();
    const double ylow = dogParamsCart3.get_ylow();
    const double zlow = dogParamsCart3.get_zlow();

    // ----------------------------------
    // Loop over all elements of interest
    // ----------------------------------    

    const int mpoints = 1;
#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    for (int j=jstart; j<=jend; j++)
    for (int k=kstart; k<=kend; k++)
    {
        const double xc = xlow + (double(i)-0.5)*dx;
        const double yc = ylow + (double(j)-0.5)*dy;
        const double zc = zlow + (double(k)-0.5)*dz;

        // each of these three items needs to be private to each thread ..
        dTensor2 xpts(mpoints,3);
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
            xpts.set( m, 3, zc );

            // Solution values (q) at each grid point
            for (int me=1; me<=meqn; me++)
            {

                for (int k=1; k<=mpoints; k++)
                {
                    qvals.set(m, me, qvals.get(m,me) + qin.get(i,j,k, me) );
                }
            }

            // Auxiliary values (aux) at each grid point
            for (int ma=1; ma<=maux; ma++)
            {
                auxvals.set(m,ma, 0.0 );
                for (int k=1; k<=mpoints; k++)
                {
                    auxvals.set(m,ma, auxvals.get(m,ma) + auxin.get(i,j,k,ma) );
                }
            }
        }

        // Call user-supplied function to set fvals
        Func(xpts, qvals, auxvals, fvals);

        // Evaluate "integrals"
        for (int m1=1; m1<=mlength; m1++)
        {
            Fout.set(i,j,k, m1, fvals.get(1,  m1) );
        }

    }

}
