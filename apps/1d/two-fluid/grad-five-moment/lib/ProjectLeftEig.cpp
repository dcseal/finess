#include <cmath>
#include "Components.h"
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that projects
// Q onto the left eigenvectors of the flux 
// Jacobian; the result is stored in W
//
class dTensor1;
class dTensor2;
void ProjectLeftEig( const dTensor1& Aux_ave,
//  const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{    

//  void ProjectLeftEig_FiveMoment( int n_offset,
//      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
//  ProjectLeftEig_FiveMoment(0, Q_ave, Q, W);
//  ProjectLeftEig_FiveMoment(5, Q_ave, Q, W);

//  void ProjectLeftEig_Maxwell(int n_offset,
//      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
//  ProjectLeftEig_Maxwell(10, Q_ave, Q, W);

//  TODO - I'm not sure what this part is for. -DS
//  void ProjectLeftConvectedScalar(int idx,
//      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
//  ProjectLeftConvectedScalar(_entropy_i, Q_ave, Q, W);
//  ProjectLeftConvectedScalar(_entropy_e, Q_ave, Q, W);

    // -- Below is copied from DoGPack -- //

    const int meqn        = Q.getsize(1);
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

    // Project onto left eigenvectors
    for(int k=1; k<=W.getsize(2); k++)
    {

        W.set(1,k, ((umag2/2.0-H-u1*c)*Q.get(2,k) 
                    + (umag2/2.0*(c-u1)+H*u1)*Q.get(1,k) 
                    + c*(Q.get(5,k)-u2*Q.get(3,k)
                        - u3*Q.get(4,k)))/(c*(2.0*H-umag2)) );

        W.set(2,k, 2.0*((H-umag2)*Q.get(1,k) + u1*Q.get(2,k) 
                    + u2*Q.get(3,k) + u3*Q.get(4,k) 
                    - Q.get(5,k))/(2.0*H-umag2) );

        W.set(3,k, Q.get(3,k)-u2*Q.get(1,k) );

        W.set(4,k, Q.get(4,k)-u3*Q.get(1,k) );

        W.set(5,k, ((H-umag2/2.0-u1*c)*Q.get(2,k) 
                    + (umag2/2.0*(c+u1)-H*u1)*Q.get(1,k)
                    + c*(Q.get(5,k)-u2*Q.get(3,k) 
                        - u3*Q.get(4,k)))/(c*(2.0*H-umag2)) );
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
    c      = sqrt(gamma*press/rho);
    H      = (energy+press)/rho; 

    // Project onto left eigenvectors
    for (int k=1; k<=W.getsize(2); k++)
    {
        W.set(6,k, ((umag2/2.0-H-u1*c)*Q.get(7,k) 
                    + (umag2/2.0*(c-u1)+H*u1)*Q.get(6,k) 
                    + c*(Q.get(10,k)-u2*Q.get(8,k)
                        - u3*Q.get(9,k)))/(c*(2.0*H-umag2)) );

        W.set(7,k, 2.0*((H-umag2)*Q.get(6,k) + u1*Q.get(7,k) 
                    + u2*Q.get(8,k) + u3*Q.get(9,k) 
                    - Q.get(10,k))/(2.0*H-umag2) );

        W.set(8,k, Q.get(8,k)-u2*Q.get(6,k) );

        W.set(9,k, Q.get(9,k)-u3*Q.get(6,k) );

        W.set(10,k, ((H-umag2/2.0-u1*c)*Q.get(7,k) 
                    + (umag2/2.0*(c+u1)-H*u1)*Q.get(6,k)
                    + c*(Q.get(10,k)-u2*Q.get(8,k) 
                        - u3*Q.get(9,k)))/(c*(2.0*H-umag2)) );
    }

    // -----------------------------
    // For the ELECTROMAGNETIC FIELD
    // -----------------------------

    // Project onto left eigenvectors
    for (int k=1; k<=W.getsize(2); k++)
    {
        W.set(11,k, (cs_light*Q.get(12,k) + Q.get(16,k))/(2.0*cs_light) );

        W.set(12,k, (cs_light*Q.get(13,k) - Q.get(15,k))/(2.0*cs_light) );

        W.set(13,k, Q.get(11,k) );

        W.set(14,k, Q.get(14,k) );

        W.set(15,k, (cs_light*Q.get(12,k) - Q.get(16,k))/(2.0*cs_light) );

        W.set(16,k, (cs_light*Q.get(13,k) + Q.get(15,k))/(2.0*cs_light) );
    }

}
