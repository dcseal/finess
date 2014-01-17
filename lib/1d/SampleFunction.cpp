#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "stdlib.h"
#include "DogParamsCart1.h"
using namespace std;

// All-purpose routine for sampling a function, and saving its data into a
// single tensor.
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensor2    node(mnodes,1)
//           dTensorBC2 auxin(1-mbc:mnodes+mbc, maux    )
//           dTensorBC2   qin(1-mbc:mnodes+mbc, meqn    )
//           dTensorBC2  Fout(1-mbc:mnodes+mbc, mlength )
// ---------------------------------------------------------------------
//
// The reason there is an extra awkward parameter, mpoints in here is to keep
// the same user interface that DoGPack uses.
//
void SampleFunction( 
    int istart, int iend,
    const dTensor2& node,
    const dTensorBC2& qin, 
    const dTensorBC2& auxin,  
    dTensorBC2& Fout,
    void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&))
{    

    const int meqn    = qin.getsize(2);
    const int maux    = auxin.getsize(2);
    const int mlength = Fout.getsize(2);

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
    const double dx   = dogParamsCart1.get_dx();
    const double xlow = dogParamsCart1.get_xlow();

    // ----------------------------------
    // Loop over all elements of interest
    // ----------------------------------    

    const int mpoints = 1;
#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    {
        double xc = xlow + (double(i)-0.5)*dx;

        // each of these three items needs to be private to each thread ..
        dTensor1 xpts(mpoints);
        dTensor2 qvals(mpoints, meqn), auxvals(mpoints, maux);
        dTensor2 fvals(mpoints, mlength);

        qvals.setall(0.);
        auxvals.setall(0.);

        // Loop over each quadrature point
        for (int m=1; m<= mpoints; m++)
        {

            // grid point x
            xpts.set( m, xc );

            // Solution values (q) at each grid point
            for (int me=1; me<=meqn; me++)
            {

                for (int k=1; k<=mpoints; k++)
                {
                    qvals.set(m,me, qvals.get(m,me) + qin.get(i,me) );
                }
            }

            // Auxiliary values (aux) at each grid point
            for (int ma=1; ma<=maux; ma++)
            {
                auxvals.set(m,ma, 0.0 );

                for (int k=1; k<=mpoints; k++)
                {
                    auxvals.set(m,ma, auxvals.get(m,ma) + auxin.get(i,ma) );
                }
            }
        }

        // Call user-supplied function to set fvals
        Func(xpts, qvals, auxvals, fvals);

        // Evaluate integrals
        for (int m1=1; m1<=mlength; m1++)
        {
            Fout.set(i, m1, fvals.get(1,  m1) );
        }

    }

}
