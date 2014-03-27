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

void Output(
        const dTensorBC2& aux,
        const dTensorBC2& q,
        double t,
        int nframe,
        string outputdir)
{

    const int mx    = q.getsize(1);
    const int meqn  = q.getsize(2);
    const int maux  = aux.getsize(2);
    const int mbc   = q.getmbc();

    // Open file -- q
    ostringstream fname1;
    fname1 << outputdir << "/" << "q" << setfill('0') 
        << setw(4) << nframe << ".dat";
    ofstream q_file(fname1.str().c_str(), ios::out );

    q_file << setprecision(16);
    q_file << setw(24) << scientific << t << endl;

    // Compute some derivatives
    void ConstructL(
            const dTensorBC2& aux,
            const dTensorBC2& q,
            dTensorBC2& Lstar,
            dTensorBC1& smax);
    dTensorBC2 qdiff(mx,meqn,mbc);
    qdiff.copyfrom(q);

    // Output each coefficient
    if( nframe > 0 )
    {
        dTensorBC1  smax(mx,mbc);
        ConstructL( aux, q, qdiff, smax );
    }

    for (int m=1; m<=meqn; m++)
    for (int i=1; i<=mx; i++)      
    {
        q_file << setw(24) << scientific << qdiff.get(i,m) << endl;
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
    for (int i=1; i<=mx; i++)      
    {
        aux_file << setw(24) << scientific << aux.get(i,m) << endl;
    }
    aux_file.close();

    // Output additional information if needed
    void Output_Extra(
            const dTensorBC2& aux,
            const dTensorBC2& q,
            double t,
            int nframe,
            string outputdir);
    Output_Extra(aux,q,t,nframe,outputdir);

}
