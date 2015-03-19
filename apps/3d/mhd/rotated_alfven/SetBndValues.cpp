#include "IniParams.h"
#include "dogdefs.h"
#include "StateVars.h"


inline void ijk_to_xyz(int i, int j, int k, double &x, double &y, double &z){
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();
    const double zlow = global_ini_params.get_zlow();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();
    const double dz = global_ini_params.get_dz();
    x = xlow + (i - 0.5) * dx;
    y = ylow + (j - 0.5) * dy;
    z = zlow + (k - 0.5) * dz;
}

inline double a1_original_to_periodic(int i, int j, int k, double a1_original){
    using std::sin;
    using std::cos;
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};
    const double b1av = n[0], b2av = n[1], b3av = n[2];
    double x, y, z;
    ijk_to_xyz(i, j, k, x, y, z);
    return a1_original - b2av * z;
}

inline double a2_original_to_periodic(int i, int j, int k, double a2_original){
    using std::sin;
    using std::cos;
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};
    const double b1av = n[0], b2av = n[1], b3av = n[2];
    double x, y, z;
    ijk_to_xyz(i, j, k, x, y, z);
    return a2_original - b3av * x;
}

inline double a3_original_to_periodic(int i, int j, int k, double a3_original){
    using std::sin;
    using std::cos;
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};
    const double b1av = n[0], b2av = n[1], b3av = n[2];
    double x, y, z;
    ijk_to_xyz(i, j, k, x, y, z);
    return a3_original - b1av * y;
}

inline double a1_periodic_to_original(int i, int j, int k, double a1_periodic){
    using std::sin;
    using std::cos;
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};
    const double b1av = n[0], b2av = n[1], b3av = n[2];
    double x, y, z;
    ijk_to_xyz(i, j, k, x, y, z);
    return a1_periodic + b2av * z;
}

inline double a2_periodic_to_original(int i, int j, int k, double a2_periodic){
    using std::sin;
    using std::cos;
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};
    const double b1av = n[0], b2av = n[1], b3av = n[2];
    double x, y, z;
    ijk_to_xyz(i, j, k, x, y, z);
    return a2_periodic + b3av * x;
}

inline double a3_periodic_to_original(int i, int j, int k, double a3_periodic){
    using std::sin;
    using std::cos;
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};
    const double b1av = n[0], b2av = n[1], b3av = n[2];
    double x, y, z;
    ijk_to_xyz(i, j, k, x, y, z);
    return a3_periodic + b1av * y;
}


// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH-ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void SetBndValues(StateVars& Q)
{

    dTensorBC4& q   = Q.ref_q();
    dTensorBC4& aux = Q.ref_aux();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int mz   = q.getsize(3);
    const int meqn = q.getsize(4);
    const int maux = aux.getsize(4);
    const int mbc  = q.getmbc();
    const double xlow = global_ini_params.get_xlow();
    const double ylow = global_ini_params.get_ylow();
    const double zlow = global_ini_params.get_zlow();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();
    const double dz = global_ini_params.get_dz();
 
    const double theta = global_ini_params.get_theta();
    const double phi = global_ini_params.get_phi();
    const double n[] = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};

    const double b1av = n[0], b2av = n[1], b3av = n[2];
   

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(0); i++)
        for (int j=1; j<=(my); j++)
            for (int k=1; k<=(mz); k++)
            {
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(mx+i,j,k,m);
                    q.set(i,j,k,m, tmp );
                }

                aux.set(i, j, k, 1,
                        a1_periodic_to_original(i, j, k,
                            a1_original_to_periodic(mx+i, j, k,
                                aux.get(mx+i, j, k, 1))));
                aux.set(i, j, k, 2,
                        a2_periodic_to_original(i, j, k,
                            a2_original_to_periodic(mx+i, j, k,
                                aux.get(mx+i, j, k, 2))));
                aux.set(i, j, k, 3,
                        a3_periodic_to_original(i, j, k,
                            a3_original_to_periodic(mx+i, j, k,
                                aux.get(mx+i, j, k, 3))));
            }
    // ***********************************************


  // ***********************************************
  // RIGHT BOUNDARY
  // ***********************************************
#pragma omp parallel for
    for (int i=(mx+1); i<=(mx+mbc); i++)
        for (int j=1; j<=(my); j++)
            for (int k=1; k<=(mz); k++)    
            {
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i-mx,j,k,m);
                    q.set(i,j,k,m, tmp );
                }
                aux.set(i, j, k, 1,
                        a1_periodic_to_original(i, j, k,
                            a1_original_to_periodic(i-mx, j, k,
                                aux.get(i-mx, j, k, 1))));
                aux.set(i, j, k, 2,
                        a2_periodic_to_original(i, j, k,
                            a2_original_to_periodic(i-mx, j, k,
                                aux.get(i-mx, j, k, 2))));
                aux.set(i, j, k, 3,
                        a3_periodic_to_original(i, j, k,
                            a3_original_to_periodic(i-mx, j, k,
                                aux.get(i-mx, j, k, 3))));
            }
    // ***********************************************


    // ***********************************************
    // FRONT BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
        for (int j=(1-mbc); j<=(0); j++)    
            for (int k=1; k<=(mz); k++)
            {

                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i,my+j,k,m);
                    q.set(i,j,k,m, tmp );
                }
                aux.set(i, j, k, 1,
                        a1_periodic_to_original(i, j, k,
                            a1_original_to_periodic(i, my+j, k,
                                aux.get(i, my+j, k, 1))));
                aux.set(i, j, k, 2,
                        a2_periodic_to_original(i, j, k,
                            a2_original_to_periodic(i, my+j, k,
                                aux.get(i, my+j, k, 2))));
                aux.set(i, j, k, 3,
                        a3_periodic_to_original(i, j, k,
                            a3_original_to_periodic(i, my+j, k,
                                aux.get(i, my+j, k, 3))));
            }
    // ***********************************************


    // ***********************************************
    // BACK BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
        for (int j=(my); j<=(my+mbc); j++)    
            for (int k=1; k<=(mz); k++)
            {
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i,j-my,k,m);
                    q.set(i,j,k,m, tmp );
                }	
                aux.set(i, j, k, 1,
                        a1_periodic_to_original(i, j, k,
                            a1_original_to_periodic(i, j-my, k,
                                aux.get(i, j-my, k, 1))));
                aux.set(i, j, k, 2,
                        a2_periodic_to_original(i, j, k,
                            a2_original_to_periodic(i, j-my, k,
                                aux.get(i, j-my, k, 2))));
                aux.set(i, j, k, 3,
                        a3_periodic_to_original(i, j, k,
                            a3_original_to_periodic(i, j-my, k,
                                aux.get(i, j-my, k, 3))));
            }
    // ***********************************************


    // ***********************************************
    // BOTTOM BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
        for (int j=(1-mbc); j<=(my+mbc); j++)    
            for (int k=(1-mbc); k<=(0); k++)
            {
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i,j,mz+k,m);
                    q.set(i,j,k,m, tmp );
                }
                aux.set(i, j, k, 1,
                        a1_periodic_to_original(i, j, k,
                            a1_original_to_periodic(i, j, mz+k,
                                aux.get(i, j, mz+k, 1))));
                aux.set(i, j, k, 2,
                        a2_periodic_to_original(i, j, k,
                            a2_original_to_periodic(i, j, mz+k,
                                aux.get(i, j, mz+k, 2))));
                aux.set(i, j, k, 3,
                        a3_periodic_to_original(i, j, k,
                            a3_original_to_periodic(i, j, mz+k,
                                aux.get(i, j, mz+k, 3))));    
            }
    // ***********************************************
  

    // ***********************************************
    // TOP BOUNDARY
    // ***********************************************
#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
        for (int j=(1-mbc); j<=(my+mbc); j++)    
            for (int k=(mz+1); k<=(mz+mbc); k++)
            {
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i,j,k-mz,m);
                    q.set(i,j,k,m, tmp );
                }
                aux.set(i, j, k, 1,
                        a1_periodic_to_original(i, j, k,
                            a1_original_to_periodic(i, j, k-mz,
                                aux.get(i, j, k-mz, 1))));
                aux.set(i, j, k, 2,
                        a2_periodic_to_original(i, j, k,
                            a2_original_to_periodic(i, j, k-mz,
                                aux.get(i, j, k-mz, 2))));
                aux.set(i, j, k, 3,
                        a3_periodic_to_original(i, j, k,
                            a3_original_to_periodic(i, j, k-mz,
                                aux.get(i, j, k-mz, 3))));    
            }
    // ***********************************************
}
