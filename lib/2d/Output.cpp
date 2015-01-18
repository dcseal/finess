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
#include "silo_wrapper.h"


namespace {
    void OutputASCII( const StateVars& Q, int nframe );
    void OutputSilo( const StateVars& Q, int nframe );
}

using namespace std;
using namespace finess;

void Output(const StateVars& Q, int nframe){
    IniParams::DataFormat::enum_type datafmt = global_ini_params.get_datafmt();
    if(datafmt == IniParams::DataFormat::ASCII)
        OutputASCII(Q, nframe);
    else if(datafmt == IniParams::DataFormat::Silo)
        OutputSilo(Q, nframe);
}


namespace {
    void OutputASCII( const StateVars& Q, int nframe )
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

    void OutputSilo( const StateVars& Q, int nframe )
    {
        const dTensorBC3& q   = Q.const_ref_q  ();
        const dTensorBC3& aux = Q.const_ref_aux();
        const double t        = Q.get_t();

        const string outputdir = global_ini_params.get_output_dir();

        const int meqn    = global_ini_params.get_meqn();
        const int maux    = global_ini_params.get_maux();
        const int mx      = global_ini_params.get_mx();
        const int my      = global_ini_params.get_my();
        const double xlow = global_ini_params.get_xlow();
        const double ylow = global_ini_params.get_ylow();
        const double dx = global_ini_params.get_dx();
        const double dy = global_ini_params.get_dy();


        ostringstream filename_oss;
        filename_oss << outputdir << "/" << "qa" << setfill('0') 
            << setw(4) << nframe << ".silo";
        const string filename = filename_oss.str();
        const string comment = "some comments.";

        silo::SetCompression setcompression_guard("METHOD=GZIP");
        silo::File qa_file(filename, comment);

        silo::TimeOptions time_opt(nframe, t);
        string meshname = "quadmesh";
        vector<double> xcoords(mx), ycoords(my);
        for(int i = 0; i < mx; ++i)
            xcoords[i] = xlow + (double(i) + 0.5) * dx;
        for(int i = 0; i < my; ++i)
            ycoords[i] = ylow + (double(i) + 0.5) * dy;
        int dims[] = {mx, my};
        int ndims = 2;
        double *coords[] = {&xcoords.front(), &ycoords.front()};
        if(DBPutQuadmesh(qa_file.get_ptr(), meshname.c_str(), NULL,
                    coords, dims, ndims, DB_DOUBLE, DB_COLLINEAR, time_opt.get_ptr())
                == -1)
            throw silo::Exception("Could not put quadmesh: nframe = " + anyToString(nframe)
                    + ", t = " + anyToString(t));

        const int nnodes = mx * my;        

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
            for(int j = 1; j <= my; ++j)
                for(int i = 1; i <= mx; ++i)
                    for(int m = 0; m < nvars; ++m)
                        v[m].push_back(q.get(i, j, m + 1));
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
            for(int j = 1; j <= my; ++j)
                for(int i = 1; i <= mx; ++i)
                    for(int m = 0; m < nvars; ++m)
                        v[m].push_back(aux.get(i, j, m + 1));
            for(int i = 0; i < nvars; ++i)
                if(DBPutQuadvar1(qa_file.get_ptr(), component_names_c[i], meshname.c_str(), components[i],
                            dims, ndims, NULL, 0,
                            DB_DOUBLE, DB_NODECENT, NULL)
                        == -1)
                    throw silo::Exception("Could not write aux: nframe = " + anyToString(nframe)
                            + ", t = " + anyToString(t));
        }
    }

}
