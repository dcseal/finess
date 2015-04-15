#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "IniParams.h"
#include "IniParams.h"
#include "dog_math.h"
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
    const double sqdx = sqrt(dx);

    string outputdir = global_ini_params.get_output_dir();
    string fname1 = outputdir+"/conservation.dat";
    string fname2 = outputdir+"/total-variation.dat";

    ofstream write_file1,write_file2;
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
    dTensor2 qsum(meqn, 4); qsum.setall(0.);
    if( global_ini_params.get_mcapa() < 1 ) // without capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            double maxq = q.get(1,m);
            for (int i=1; i<=mx; i++)
            {
                const double x     = xlow + (double(i)-0.5)*dx;
                const double qtmp  = q.get(i, m);
                const double abs_q = fabs( q.get(i, m) );

                qsum.set(m, 1, qsum.get(m,1) +   dx*qtmp       );  // total mass
                qsum.set(m, 2, qsum.get(m,2) +   dx*abs_q      );  // L1-norm
                qsum.set(m, 3, qsum.get(m,3) + sqdx*qtmp*qtmp  );  // L2-norm
                qsum.set(m, 4, Max( abs_q, qsum.get(m,4) )     );  // L-inf norm

            }
            qsum.set( m, 3, sqrt( qsum.get(m,3) ) );
        }
    }
    else // with capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            for (int i=1; i<=mx; i++)
            {
                const double x    = xlow + (double(i)-0.5)*dx;
                const double qtmp = q.get(i,m);
                const double atmp = aux.get(i, global_ini_params.get_mcapa() );
                const double abs_q = fabs( atmp*q.get(i, m) );

                qsum.set(m, 1, qsum.get(m,1) + atmp*dx*qtmp        );  // total mass
                qsum.set(m, 2, qsum.get(m,2) + atmp*dx*abs_q       );  // L1-norm
                qsum.set(m, 3, qsum.get(m,3) + atmp*sqdx*qtmp*qtmp );  // L2-norm
                qsum.set(m, 4, Max( abs_q, qsum.get(m,4) )         );  // L-inf norm
            }
            qsum.set( m, 3, sqrt( qsum.get(m,3) ) );
        }
    }

    // Write data to file
    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        write_file1 << setw(24) << scientific;
        for( int s=1; s <= 4; s++ )
        {
            write_file1 << qsum.get(m,s) << "  ";
        }
    }
    write_file1 << endl;


    // Total variation
    dTensor1 tv(meqn);
    for (int m=1; m<=meqn; m++)
    {
        tv.set(m, fabs(q.get(1,m) - q.get(mx,m) ) );
        for (int i=2; i<=mx; i++)
        {
            tv.set(m, tv.get(m) + fabs( q.get(i,m) - q.get(i-1,m) ) );
        }
    }

    // Maximum and minimum value observed during simulation
    dTensor1 MaxQ(meqn); 
    dTensor1 MinQ(meqn);
    for (int m=1; m<=meqn; m++)
    {
        MaxQ.set(m, q.get(1,m) );
        MinQ.set(m, q.get(1,m) );
        for (int i=2; i<=mx; i++)
        {
            MaxQ.set(m, Max( MaxQ.get(m), q.get(i,m) ) );
            MinQ.set(m, Min( MinQ.get(m), q.get(i,m) ) );
        }
    }


    write_file2 << setprecision(16);
    write_file2 << setw(24) << scientific << t << " ";
    for (int m=1; m<=meqn; m++)
    {
        write_file2 << setw(24) << scientific << tv.get(m) << " ";
        write_file2 << setw(24) << scientific << MinQ.get(m) << " ";
        write_file2 << setw(24) << scientific << MaxQ.get(m) << " ";
    }
    write_file2 << endl;

}
