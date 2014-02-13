///@file lib/1d/Output.cpp
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

/**
@brief Outputs snapshot to file.

@param node Not used.
@param aux Array to be output to <tt><outputdir>/a\<nframe\>.dat</tt>
@param q Array to be output to <tt><outputdir>/q\<nframe\>.dat</tt>
@param t Supposed to be time
@param nframe <tt>\<nframe\></tt> := <tt>nframe</tt> padded with zero
@param outputdir Path to output files

Does two things:
- Outputs to <tt><outputdir>/q\<nframe\>.dat</tt>
-- <tt>t</tt>
-- <tt>q</tt>, with 2nd index being the row index, and 1st being column
- Outputs to <tt><outputdir>/a\<nframe\>.dat</tt>
-- <tt>t</tt>
-- <tt>aux</tt>, with 2nd idnex being the row index, and 1st being column

@sa Output_Extra(...)

@todo Clarify the intent of reserving <tt>node</tt> in the arguments.
*/
void Output(const dTensor2& node, 
        const dTensorBC2& aux,
        const dTensorBC2& q,
        double t,
        int nframe,
        string outputdir)
{

    const int melems  = q.getsize(1);
    const int meqn    = q.getsize(2);
    const int maux    = aux.getsize(2);

    // Open file -- q
    ostringstream fname1;
    fname1 << outputdir << "/" << "q" << setfill('0') 
        << setw(4) << nframe << ".dat";
    ofstream q_file(fname1.str().c_str(), ios::out );

    q_file << setprecision(16);
    q_file << setw(24) << scientific << t << endl;

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

    // Output additional information if needed - TODO reintroduce this call
    void Output_Extra(const dTensor2& node, 
            const dTensorBC2& aux,
            const dTensorBC2& q,
            double t,
            int nframe,
            string outputdir);
    Output_Extra(node,aux,q,t,nframe,outputdir);

}
