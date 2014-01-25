#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "DogParams.h"

using namespace std;

void ConSoln( const dTensor2& node, 
    const dTensorBC2& aux,
    const dTensorBC2& q, 
    double t, string outputdir)
{

    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);
    string fname1 = outputdir+"/conservation.dat";
    ofstream write_file1,write_file2;
    dTensor1 qsum(meqn);
    dTensor1 res_sum(meqn);

    if (t==0) 
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

            for (int i=1; i<=melems; i++)
            {
                double x = node.get(i,1);
                double dtmp = node.get(i+1,1)-node.get(i,1);
                double qtmp = q.get(i, m);

                qsum.set(m, (qsum.get(m) + dtmp*qtmp) );
            }
        }
    }
    else // with capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m, 0.0);

            for (int i=1; i<=melems; i++)
            {
                double x = node.get(i,1);
                double dtmp = node.get(i+1,1)-node.get(i,1);
                double qtmp = q.get(i,m);
                double atmp = aux.get(i, dogParams.get_mcapa() );

                qsum.set(m, (qsum.get(m) + atmp*dtmp*qtmp) );
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
