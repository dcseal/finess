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

    const dTensorBC2& q   = Q.const_ref_q  ();
    const dTensorBC2& aux = Q.const_ref_aux();
    const double t        = Q.get_t();

    const int melems  = q.getsize(1);
    const int meqn    = q.getsize(2);
    const int maux    = aux.getsize(2);

    string outputdir = global_ini_params.get_output_dir();

    // Open file -- q
    ostringstream fname1;
    fname1 << outputdir << "/" << "q" << setfill('0') 
        << setw(4) << nframe << ".dat";
    ofstream q_file(fname1.str().c_str(), ios::out );

    q_file << setprecision(16);
    q_file << setw(24) << scientific << Qstate.get_t() << endl;

    // Output each coefficient
    for (int m=1; m<=meqn; m++)
    for (int i=1; i<=melems; i++)      
    {
        q_file << setw(24) << scientific << q.get(i,m) << endl;
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
    for (int i=1; i<=melems; i++)      
    {
        aux_file << setw(24) << scientific << aux.get(i,m) << endl;
    }
    aux_file.close();

    // Output additional information (if relinked)
    void Output_Extra( const StateVars& Qstate, int nframe );
    Output_Extra(Qstate, nframe );

}
