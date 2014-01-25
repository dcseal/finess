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

    const int     mx = dogParamsCart2.get_mx();
    const int     my = dogParamsCart2.get_my();
    const int   meqn = dogParams.get_meqn();
    const int   maux = dogParams.get_maux();

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
//          qsum.set(m,0.0);

//          for (int i=1; i<=mx; i++)
//          for (int j=1; j<=my; j++)
//          {
//              double x = node.get(i,j,1);
//              double dtmp = node.get(i+1,1)-node.get(i,1);
//              double qtmp = q.get(i, m);

// TODO
                qsum.set(m, 0.0 );
//          }
        }
    }
    else // with capacity function
    {
        for (int m=1; m<=meqn; m++)
        {
            qsum.set(m, 0.0);

//          for (int i=1; i<=mx; i++)
//          for (int j=1; j<=my; j++)
//          {
//              double x = node.get(i,j,1);
//              double y = node.get(i,j,2);
//              double dtmp = node.get(i+1,1)-node.get(i,1);

//              double qtmp = q.get(i,j,m);
//              double atmp = aux.get(i,j, dogParams.get_mcapa() );

//              qsum.set(m, (qsum.get(m) + atmp*dtmp*qtmp) );
// TODO
//          }
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
