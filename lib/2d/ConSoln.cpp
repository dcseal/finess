#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"

using namespace std;

void ConSoln( const dTensorBC3& aux, const dTensorBC3& q, 
              double t, string outputdir)
{

    // Size of the solution
    const int     mx = dogParamsCart2.get_mx();
    const int     my = dogParamsCart2.get_my();
    const int   meqn = dogParams.get_meqn();
    const int   maux = dogParams.get_maux();

    // Grid information:
    const double dx   = dogParamsCart2.get_dx();
    const double dy   = dogParamsCart2.get_dy();
    const double xlow = dogParamsCart2.get_xlow();
    const double ylow = dogParamsCart2.get_ylow();

    string fname1 = outputdir+"/conservation.dat";
    ofstream write_file1,write_file2;
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
            for (int j=1; j<=my; j++)
            {
                const double x    = xlow + (double(i)-0.5)*dx;
                const double y    = ylow + (double(j)-0.5)*dy;
                const double qtmp = q.get(i,j,m);
                qsum.set(m, qsum.get(m) + dx*dy*qtmp );
            }
        }
    }
    else // with capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m, 0.0);

            for (int i=1; i<=mx; i++)
            for (int j=1; j<=my; j++)
            {
                const double x  = xlow + (double(i)-0.5)*dx;
                const double y  = ylow + (double(j)-0.5)*dy;

                double qtmp = q.get(i,j,m);
                double atmp = aux.get(i,j, dogParams.get_mcapa() );
                qsum.set(m, (qsum.get(m) + atmp*dx*dy*qtmp) );
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
