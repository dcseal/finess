#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"

using namespace std;

void ConSoln( const StateVars& Qstate )
{

    const dTensorBC2& q   = Qstate.const_ref_q  ();
    const dTensorBC2& aux = Qstate.const_ref_aux();
    const double        t = Qstate.get_t();

    const int     mx = global_ini_params.get_mx();
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();

    // Grid information:
    const double dx   = global_ini_params.get_dx();
    const double xlow = global_ini_params.get_xlow();

    string outputdir = global_ini_params.get_output_dir();
    string fname1 = outputdir+"/conservation.dat";
    string fname2 = outputdir+"/total-variation.dat";

    ofstream write_file1,write_file2;
    dTensor1 qsum(meqn);
    dTensor1 res_sum(meqn);

    if( t==0 ) 
    {
        write_file1.open(fname1.c_str(), ofstream::out);
        write_file2.open(fname2.c_str(), ofstream::out);
    }
    else
    {
        write_file1.open(fname1.c_str(), ofstream::app);
        write_file2.open(fname2.c_str(), ofstream::app);
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
            {
                const double x    = xlow + (double(i)-0.5)*dx;
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
                const double atmp = aux.get(i, global_ini_params.get_mcapa() );

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


    // Total variation
    for (int m=1; m<=meqn; m++)
    {
        qsum.set(m, fabs(q.get(1,m) - q.get(mx,m) ) );
        for (int i=2; i<=mx; i++)
        {
            qsum.set(m, qsum.get(m) + fabs( q.get(i,m) - q.get(i-1,m) ) );
//          printf("TF(%d) = %f\n", i, fabs(q.get(i,m)-q.get(i-1,m) ) );
        }
    }

    write_file2 << setprecision(16);
    write_file2 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        write_file2 << setw(24) << scientific << qsum.get(m) << " ";
    }
    write_file2 << endl;

}
