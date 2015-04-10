#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"

#ifdef USE_SILO
#include "silo_wrapper.h"
#endif

namespace {
    void OutputASCII( const StateVars& Q, int nframe );
#ifdef USE_SILO    
    void OutputSilo( const StateVars& Q, int nframe );
#endif
}

using namespace std;
#ifdef USE_SILO
using namespace finess;
#endif

void Output(const StateVars& Q, int nframe){
    IniParams::DataFormat::enum_type datafmt = global_ini_params.get_datafmt();
    if(datafmt == IniParams::DataFormat::ASCII)
        OutputASCII(Q, nframe);
    else if(datafmt == IniParams::DataFormat::Silo){
#ifdef USE_SILO        
        OutputSilo(Q, nframe);
#else
        cerr << "Silo format is not supported.  Using ASCII instead." << endl;
        OutputASCII(Q, nframe);
#endif
    }
}

namespace {
    void OutputASCII( const StateVars& Q, int nframe )
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
        q_file << setw(24) << scientific << Q.get_t() << endl;

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
        //        void Output_Extra( const StateVars& Q, int nframe );
        //        Output_Extra(Q, nframe );
    }
#ifdef USE_SILO
    void OutputSilo( const StateVars& Q, int nframe )
    {
        const dTensorBC2& q   = Q.const_ref_q  ();
        const dTensorBC2& aux = Q.const_ref_aux();
        const double t        = Q.get_t();

        const string outputdir = global_ini_params.get_output_dir();

        const int meqn    = global_ini_params.get_meqn();
        const int maux    = global_ini_params.get_maux();
        const int mx      = global_ini_params.get_mx();
        const double xlow = global_ini_params.get_xlow();
        const double dx = global_ini_params.get_dx();


        ostringstream filename_oss;
        filename_oss << outputdir << "/" << "qa" << setfill('0') 
            << setw(4) << nframe << ".silo";
        const string filename = filename_oss.str();
        const string comment = "some comments.";

        silo::SetCompression setcompression_guard("METHOD=GZIP");
        silo::File qa_file(filename, comment);

        silo::TimeOptions time_opt(nframe, t);
        string meshname = "quadmesh";
        vector<double> xcoords(mx);
        for(int i = 0; i < mx; ++i)
            xcoords[i] = xlow + (double(i) + 0.5) * dx;
        int dims[] = {mx};
        int ndims = 1;
        double *coords[] = {&xcoords.front()};
        if(DBPutQuadmesh(qa_file.get_ptr(), meshname.c_str(), NULL,
                    coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, time_opt.get_ptr())
                == -1)
            throw silo::Exception("Could not put quadmesh: nframe = " + anyToString(nframe)
                    + ", t = " + anyToString(t));

        const int nnodes = mx;        

        {
            string varname = "q";
            int nvars = meqn;
            string component_names[nvars];
            char* component_names_c[nvars];
            for(int i = 0; i < nvars; ++i){
                component_names[i] = string("q") + anyToString(i + 1);
                component_names_c[i] = const_cast<char*>(component_names[i].c_str());
            }
            double* components[nvars];
            vector<double> v[nvars];
            for(int i = 0; i < nvars; ++i){
                v[i].reserve(nnodes);
                components[i] = &v[i].front();
            }
            for(int i = 1; i <= mx; ++i)
                for(int m = 0; m < nvars; ++m)
                    v[m].push_back(q.get(i, m + 1));

            for(int i = 0; i < nvars; ++i)
                if(DBPutQuadvar1(qa_file.get_ptr(), component_names_c[i], meshname.c_str(), components[i], 
                            dims, ndims, NULL, 0,
                            DB_DOUBLE, DB_NODECENT, NULL)
                        == -1)
                    throw silo::Exception("Could not write q: nframe = " + anyToString(nframe)
                            + ", t = " + anyToString(t));

        }

        if(maux > 0){
            string varname = "a";
            int nvars = maux;
            string component_names[nvars];
            char* component_names_c[nvars];
            for(int i = 0; i < nvars; ++i){
                component_names[i] = string("a") + anyToString(i + 1);
                component_names_c[i] = const_cast<char*>(component_names[i].c_str());
            }
            double* components[nvars];
            vector<double> v[nvars];
            for(int i = 0; i < nvars; ++i){
                v[i].reserve(nnodes);
                components[i] = &v[i].front();
            }
            for(int i = 1; i <= mx; ++i)
                for(int m = 0; m < nvars; ++m)
                    v[m].push_back(aux.get(i, m + 1));
            for(int i = 0; i < nvars; ++i)
                if(DBPutQuadvar1(qa_file.get_ptr(), component_names_c[i], meshname.c_str(), components[i], 
                            dims, ndims, NULL, 0,
                            DB_DOUBLE, DB_NODECENT, NULL)
                        == -1)
                    throw silo::Exception("Could not write aux: nframe = " + anyToString(nframe)
                            + ", t = " + anyToString(t));
        }
    }
#endif
}
