///@file lib/1d/ConSoln.cpp

#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart1.h"

using namespace std;

///@brief Outputs the conserved quantities to file.
///
///Appends (unless <tt>t</tt>=0, in which case creates and writes) to <tt>\<outputdir\>/conservation.dat</tt> a line: <br>
/// \<t\>             @f$ (\sum_q^{\mathrm{mx}}(q_m)_i C_i \mathrm{dx})_{m=1, \ldots, \mathrm{meqn}} @f$  <br>
///where @f$ q_m @f$ is the m-th column of <tt>q</tt>,
///and @f$ C_i @f$ is the capacity function values, which are 1 if dogParams.get_mcapa() returns 0, and <tt>aux[i, dogParams.get_mcapa()]</tt> otherwise.
void ConSoln( 
    const dTensorBC2& aux,
    const dTensorBC2& q, 
    double t, string outputdir)
{

    const int     mx = dogParamsCart1.get_mx();
    const int   meqn = dogParams.get_meqn();
    const int   maux = dogParams.get_maux();

    // Grid information:
    const double dx   = dogParamsCart1.get_dx();
    const double xlow = dogParamsCart1.get_xlow();

    string fname1 = outputdir+"/conservation.dat";
    ofstream write_file1,write_file2;   //write_file2 is not used!!
    dTensor1 qsum(meqn);
    dTensor1 res_sum(meqn);

    if( t==0 ) 
    {
        write_file1.open(fname1.c_str(), ofstream::out);
    }
    else
    {
        write_file1.open(fname1.c_str(), ofstream::app);
    }

    // -----------------
    // CONSERVATION
    // -----------------
    if( dogParams.get_mcapa() < 1 ) // without capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m,0.0);

            for (int i=1; i<=mx; i++)
            {
                const double x    = xlow + (double(i)-0.5)*dx;   //not used!!
                const double qtmp = q.get(i, m);

                qsum.set(m, qsum.get(m) + dx*qtmp );
            }
        }
    }
    else // with capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m, 0.0);

            for (int i=1; i<=mx; i++)
            {
                const double x    = xlow + (double(i)-0.5)*dx;
                const double qtmp = q.get(i,m);
                const double atmp = aux.get(i, dogParams.get_mcapa() );

                qsum.set(m, (qsum.get(m) + atmp*dx*qtmp) );
            }
        }
    }

    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        if (abs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
        write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
    write_file1 << endl;

}
