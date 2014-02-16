#include <iostream>
#include <cmath>
#include "dogdefs.h"
#include "stdlib.h"
#include "DogParamsCart1.h"
using namespace std;


///@brief  All-purpose routine for sampling a function, and saving its data into a
///        single tensor.
///
///@param istart  
///@param iend
///@param node   An array that never gets used.
///@param qin    An array that is passed on to <tt>Func</tt>
///@param auxin  An array that is passed on to <tt>Func</tt>
///@param Fout   An array to hold the sampled function values.
///@param Func   The function to sample.
///
/// After a call to this function, <tt>Fout</tt> will be filled with the value of a function
/// (implemented in <tt>Func</tt>),
/// the input to which is
/// - <tt>xlow+(i-0.5)*dx</tt>, where <tt>i = istart, ..., iend</tt>
/// - <tt>qin</tt>
/// - <tt>auxin</tt>
/// 
/// The output of <tt>Func</tt> will be stored in <tt>Fout</tt>.
///
///
///*********************************
/// (Old comments) Inputs should have the following sizes:    <br>
///           dTensor2    node(mnodes,1)   <br>
///           dTensorBC2 auxin(1-mbc:mnodes+mbc, maux    ) <br>
///           dTensorBC2   qin(1-mbc:mnodes+mbc, meqn    ) <br>
///           dTensorBC2  Fout(1-mbc:mnodes+mbc, mlength ) <br>
/// The reason there is an extra awkward parameter, mpoints in here is to keep
/// the same user interface that DoGPack uses.
///
///@note Implementation details:
//    The funny logic here:
///   - Func(...) is supposed to be an function from (several) arrays to an array;
///   - However, Func(...) is not called like an array to array function.
///   -- The outer 'for' loop (i=istart, ..., iend) breaks the array into row vectors,
///      i.e. entries with same first index are grouped into a one-dimensional array.
///   -- Then these one-dimensional arrays are wrapped into a two-dimensional array,
///      whose first dimension has length 1.
///   -- These two-dimensional arrays of trivial first dimension length
///      is what get passed to <tt>Func</tt>.
///   -- The sink argument passed to <tt>Func</tt> also has trivial first dimension length.
///      Thus the result from <tt>Func</tt> is a one-dimensional array wrapped in a two-dimensional array.
///   -- This one-dimensional array is to be interpreted 
///      as a row vector in the two-dimensional array <tt>Fout</tt>.
///   -- Such an interpretation is materialized as the last statement of the main 'for' loop.
///      
///     
///@note Collaboration details: This function is called three times:
///- In ConstructL(...), where <tt>Func</tt> is SourceTermFunc(...),
///  a user-supplied function that takes four arguments;
///- In RunFinpack(...), where <tt>Func</tt> is AuxFunc(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&),
///  a wrapper function to a user-supplied function AuxFunc(const dTensor1&, dTensor2&)
///- In RunFinpack(...), where <tt>Func</tt> is QinitFunc(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&),
///  a wrapper function to a user-supplied function QinitFunc(const dTensor1&, dTensor2&)
///@note Implementation details:  the local variable <tt>mpoints</tt> is set to be constant 1.
///
///@todo Finish this documentation.
///@todo Clarify use of <tt>node</tt> -- it is never used in this function!!
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
