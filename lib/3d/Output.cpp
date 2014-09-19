#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "IniParams.h"
#include "IniParams.h"
using namespace std;

void Output(
        const dTensorBC4& aux,
        const dTensorBC4& q,
        double t,
        int nframe,
        string outputdir)
{

    const int meqn    = global_ini_params.get_meqn();
    const int maux    = global_ini_params.get_maux();
    const int mx      = global_ini_params.get_mx();
    const int my      = global_ini_params.get_my();
    const int mz      = global_ini_params.get_mz();

    // Open file -- q
    ostringstream fname1;
    fname1 << outputdir << "/" << "q" << setfill('0') 
        << setw(4) << nframe << ".dat";
    ofstream q_file(fname1.str().c_str(), ios::out );

    q_file << setprecision(16);
    q_file << setw(24) << scientific << t << endl;

    // Output each coefficient - TODO, we could potentially reverse the order
    // here, but then the plotting routines will have to change as well.  It
    // is faster to index the arrays in odometer order. (-DS).
    for (int m=1; m<=meqn; m++)
    for (int k=1; k<=mz; k++)      
    for (int j=1; j<=my; j++)      
    for (int i=1; i<=mx; i++)      
    {
    //TODO
        q_file << setw(24) << scientific << q.get(i,j,k,m) << endl;
    }
    q_file.close();

    // Open file -- aux
    ostringstream fname2;
    fname2 << outputdir << "/" << "a" << setfill('0') 
        << setw(4) << nframe << ".dat";
    ofstream aux_file(fname2.str().c_str(), ios::out );

    aux_file << setprecision(16);
    aux_file << setw(24) << scientific << t << endl;

    // Output aux array
    for (int m=1; m<=maux; m++)
    for (int k=1; k<=mz; k++)      
    for (int j=1; j<=my; j++)      
    for (int i=1; i<=mx; i++)      
    {
        aux_file << setw(24) << scientific << aux.get(i,j,k,m) << endl;
    }
    aux_file.close();

    // Output additional information if needed - TODO reintroduce this call
//  void Output_Extra(const dTensor2& node, 
//          const dTensorBC2& aux,
//          const dTensorBC2& q,
//          double t,
//          int nframe,
//          string outputdir);
//  Output_Extra(node,aux,q,t,nframe,outputdir);

}
