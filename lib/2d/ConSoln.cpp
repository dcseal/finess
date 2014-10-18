#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "StateVars.h"

#include "IniParams.h"
using namespace std;

void ConSoln( const StateVars& Q )
{

    const dTensorBC3&    q = Q.const_ref_q  ();
    const dTensorBC3&  aux = Q.const_ref_aux();
    const double         t = Q.get_t();

    string outputdir = global_ini_params.get_output_dir();

    // Size of the solution
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();

    // Grid information:
    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();

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
    if( global_ini_params.get_mcapa() < 1 ) // without capacity function
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
                double atmp = aux.get(i,j, global_ini_params.get_mcapa() );
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
