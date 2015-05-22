#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"
#include "Output.h"

using namespace std;

void ConSoln( const StateVars& Q )
{

    const dTensorBC4& q   = Q.const_ref_q  ();
    const dTensorBC4& aux = Q.const_ref_aux();
    const double t        = Q.get_t();
    
    string outputdir = global_ini_params.get_output_dir();

    // Size of the solution
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int     mz = global_ini_params.get_mz();
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();

    // Grid information:
    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();
    const double dz   = global_ini_params.get_dz();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();
    const double zlow = global_ini_params.get_zlow();

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
            for (int k=1; k<=mz; k++)
            {
                const double x    = xlow + (double(i)-0.5)*dx;
                const double y    = ylow + (double(j)-0.5)*dy;
                const double z    = zlow + (double(k)-0.5)*dz;

// TODO
                const double qtmp = q.get(i,j,k,m);
                qsum.set(m, qsum.get(m) + dx*dy*dz*qtmp );
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
            for (int k=1; k<=mz; k++)
            {
                const double x  = xlow + (double(i)-0.5)*dx;
                const double y  = ylow + (double(j)-0.5)*dy;
                const double z  = zlow + (double(k)-0.5)*dz;

                double qtmp = q.get(i,j,k,m);
                double atmp = aux.get(i,j,k, global_ini_params.get_mcapa() );
                qsum.set(m, (qsum.get(m) + atmp*dx*dy*dz*qtmp) );
            }
        }
    }

    write_file1 << setprecision(16);
    write_file1 << setw(24) << scientific << t << " ";
//    for (int m=1; m<=meqn; m++)
//    {
//        if (abs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
//        write_file1 << setw(24) << scientific << qsum.get(m) << " ";
//    }
    double gamma = global_ini_params.get_gamma();
    using std::min;
    double min_p = std::numeric_limits<double>::infinity();
#pragma omp parallel for reduction(min:min_p)
    for (int i=1; i<=mx; i++)
	for (int j=1; j<=my; j++)
	    for (int k=1; k<=mz; k++){
	    	double rho = q.get(i, j, k, 1);
		double u1 = q.get(i, j, k, 2) / rho;
		double u2 = q.get(i, j, k, 3) / rho;
		double u3 = q.get(i, j, k, 4) / rho;
		double E = q.get(i, j, k, 5);
		double B1 = q.get(i, j, k, 6);
		double B2 = q.get(i, j, k, 7);
		double B3 = q.get(i, j, k, 8);
		double p = (gamma - 1.0) * (E - 0.5 * rho * (u1*u1 + u2*u2 + u3*u3) - 0.5 * (B1*B1 + B2*B2 + B3*B3));
		min_p = min(min_p, p);
	    }
 
    write_file1 << setw(24) << scientific << min_p;
    
    write_file1 << endl;

    int nout = global_ini_params.get_nout();
    if(min_p < 0.0)
	Output(Q, nout + 1);


}
