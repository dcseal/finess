#include <cmath>
#include "Components.h"
#include "dogdefs.h"
#include "IniParams.h"


// This is a user-supplied routine that projects
// W onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Q
//
class dTensor1;
class dTensor2;
void ProjectRightEig( const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{    

//  void ProjectRightEig_FiveMoment(int n_offset,
//      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
//  ProjectRightEig_FiveMoment(0, Q_ave, W, Q);
//  ProjectRightEig_FiveMoment(5, Q_ave, W, Q);

//  void ProjectRightEig_Maxwell(int n_offset,
//      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
//  ProjectRightEig_Maxwell(10, Q_ave, W, Q);

//  TODO - I'm not sure what this part is for.  -DS
//  void ProjectRightConvectedScalar(int idx,
//      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
//  ProjectRightConvectedScalar(_entropy_i, Q_ave, W, Q);
//  ProjectRightConvectedScalar(_entropy_e, Q_ave, W, Q);


    // -- Below was copied straight from DoGPack -- //

    const int meqn = Q.getsize(1);
    const double cs_light = global_ini_params.get_cs_light();

    // -----------------------------
    // For the IONS
    // -----------------------------

    // Average states
    const double& gamma  = global_ini_params.get_gamma();
    double rho    = Q_ave.get(1);
    double u1     = Q_ave.get(2)/rho;
    double u2     = Q_ave.get(3)/rho;
    double u3     = Q_ave.get(4)/rho;
    double energy = Q_ave.get(5);
    double umag2  = (u1*u1 + u2*u2 + u3*u3);
    double press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
    double c      = sqrt(gamma*press/rho);
    double H      = (energy+press)/rho; 


    // Project onto right eigenvectors
    for(int k=1; k<=Q.getsize(2); k++)
    {
        Q.set(1,k, W.get(1,k) + W.get(2,k) + W.get(5,k)  );

        Q.set(2,k, (u1-c)*W.get(1,k) + u1*W.get(2,k) 
                + (u1+c)*W.get(5,k) );

        Q.set(3,k, u2*(W.get(1,k) + W.get(2,k) + W.get(5,k))
                + W.get(3,k) );

        Q.set(4,k, u3*(W.get(1,k) + W.get(2,k) + W.get(5,k))
                + W.get(4,k) );

        Q.set(5,k, (H-u1*c)*W.get(1,k) + umag2/2.0*W.get(2,k)
                + u2*W.get(3,k) + u3*W.get(4,k) 
                + (H+u1*c)*W.get(5,k) );
    }


    // -----------------------------
    // For the ELECTRONS
    // -----------------------------

    // Average states
    rho    = Q_ave.get(6);
    u1     = Q_ave.get(7)/rho;
    u2     = Q_ave.get(8)/rho;
    u3     = Q_ave.get(9)/rho;
    energy = Q_ave.get(10);
    umag2  = (u1*u1 + u2*u2 + u3*u3);
    press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
    c      = sqrt(fabs(gamma*press/rho));
    H      = (energy+press)/rho;  

    // Project onto right eigenvectors
    for(int k=1; k<=Q.getsize(2); k++)
    {
        Q.set(6,k, W.get(6,k) + W.get(7,k) + W.get(10,k)  );

        Q.set(7,k, (u1-c)*W.get(6,k) + u1*W.get(7,k) 
                + (u1+c)*W.get(10,k) );

        Q.set(8,k, u2*(W.get(6,k) + W.get(7,k) + W.get(10,k))
                + W.get(8,k) );

        Q.set(9,k, u3*(W.get(6,k) + W.get(7,k) + W.get(10,k))
                + W.get(9,k) );

        Q.set(10,k, (H-u1*c)*W.get(6,k) + umag2/2.0*W.get(7,k)
                + u2*W.get(8,k) + u3*W.get(9,k) 
                + (H+u1*c)*W.get(10,k) );
    }


    // -----------------------------
    // For the ELECTROMAGNETIC FIELD
    // -----------------------------

    // Project onto right eigenvectors
    for(int k=1; k<=Q.getsize(2); k++)
    {
        Q.set(11,k, W.get(13,k) );

        Q.set(12,k, W.get(11,k) + W.get(15,k) );

        Q.set(13,k, W.get(12,k) + W.get(16,k) );

        Q.set(14,k, W.get(14,k) );

        Q.set(15,k, cs_light*( W.get(16,k) - W.get(12,k) ) );

        Q.set(16,k, cs_light*( W.get(11,k) - W.get(15,k) ) );
    }

}
