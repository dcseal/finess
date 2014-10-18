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

void Output( const StateVars& Q, int nframe )
{

    const dTensorBC3& q   = Q.const_ref_q  ();
    const dTensorBC3& aux = Q.const_ref_aux();
    const double t        = Q.get_t();

    const int meqn    = global_ini_params.get_meqn();
    const int maux    = global_ini_params.get_maux();
    const int mx      = global_ini_params.get_mx();
    const int my      = global_ini_params.get_my();

    string outputdir = global_ini_params.get_output_dir();

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
    for (int j=1; j<=my; j++)      
    for (int i=1; i<=mx; i++)      
    {
        q_file << setw(24) << scientific << q.get(i,j,m) << endl;
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
    for (int j=1; j<=my; j++)      
    for (int i=1; i<=mx; i++)      
    {
        aux_file << setw(24) << scientific << aux.get(i,j,m) << endl;
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
