#include <cmath>
#include <iostream>
#include "stdio.h"
#include "dog_math.h"
#include "IniParams.h"
#include "constants.h"
#include "tensors.h"
using namespace std;

// MPP limiter for 2D Euler equations.
//
// This limiter works by considering a high order flux, (fHat, gHat), and
// low-order flux, (fLF, gLF), and then creates a new flux of the form:
//
//     fhat = theta fhat + (1-theta) fLF
//     ghat = theta ghat + (1-theta) gLg
//
// so that the conservative update, 
//
// q = q - dt/dx( f_{i+1/2} - f_{i-1/2} ) - dt/dy( g_{j+1/2} - g_{j-1/2} )
//
// is conservative.
//
// See: "http://arxiv.org/abs/1411.0328" and references therein for more
// details.
void ApplyMPPLimiter2D( 
        const double dt, const dTensorBC3& q, 
        const dTensorBC3& fLF, const dTensorBC3& gLF,
        dTensorBC3& fHat, dTensorBC3& gHat )
{

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();

    // Parameters used for limiting
    const double eps = 1.0e-12;
    const double dtx = dt/dx, dty = dt/dy;
    const double rho_min = eps;

    // Parameters specific for Euler equations
    const double gamma = global_ini_params.get_gamma();
    const double gm1 = gamma - 1.0;

    // Storage
    double gmin[mx][my],amin[mx][my][4],thex[mx+1][my],they[mx][my+1],ff[4],aa[4];
    double qh[meqn],qtmp[meqn],qtmp2[meqn],qtmp3[meqn],qlf[meqn],fhat_local[4][meqn],flf_local[4][meqn];
    double rho,u1,u2,u3,energy,plow,phigh,ffsum,ftmp,phigh2;
    double ac,bc,cc,delta,root1,root2,rate;
    int isgn[4];

    /////////////////////////////////////////
    // Step 1: limit rho
    /////////////////////////////////////////
    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
        gmin[i-1][j-1] = rho_min - ( q.get(i,j,1) 
                       + dtx*( fLF.get(i,j,1) - fLF.get(i+1,j,1))
                       + dty*( gLF.get(i,j,1) - gLF.get(i,j+1,1)));
 
    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
    {
        ff[0] = dtx*(fHat.get(i,j,1)  -fLF.get(i,j,1)  );
        ff[1] =-dtx*(fHat.get(i+1,j,1)-fLF.get(i+1,j,1));
        ff[2] = dty*(gHat.get(i,j,1)  -gLF.get(i,j,1)  );
        ff[3] =-dty*(gHat.get(i,j+1,1)-gLF.get(i,j+1,1));

        for (int k=0; k<4; k++)
        {
            if (ff[k]<0.0)
                isgn[k] = 1;
            else
                isgn[k] = 0;
        }

        ffsum = double(isgn[0])*ff[0];

        for (int k=1; k<4; k++)
            ffsum = ffsum + double(isgn[k])*ff[k];

        for (int k=0; k<4; k++)
        {
            if (isgn[k]==1)
                amin[i-1][j-1][k] = Min(gmin[i-1][j-1]/(ffsum-eps),1.e0);
            else
                amin[i-1][j-1][k] = 1.e0;
        }
    }

    /////////////////////////////////////////
    // Step 2: limit pressure
    /////////////////////////////////////////
    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
    {

        // Write the update for a LF solve
        for( int m=1; m <= meqn; m++ )  
            qlf[m-1]=q.get(i,j,m) + dtx*(fLF.get(i,j,m)-fLF.get(i+1,j,m))
                                  + dty*(gLF.get(i,j,m)-gLF.get(i,j+1,m));

        rho    = qlf[0];
        u1     = qlf[1]/rho;
        u2     = qlf[2]/rho;
        u3     = qlf[3]/rho;
        energy = qlf[4];
        plow   = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

        // Quick error check
        if (plow < 0.0)
            cout << "Negative solution in Lax-Friedrich's flux" << endl;

        for( int m=1; m <= meqn; m++ )
        {
            fhat_local[0][m-1] = fHat.get(i,j,m);
            fhat_local[1][m-1] = fHat.get(i+1,j,m);
            fhat_local[2][m-1] = gHat.get(i,j,m);
            fhat_local[3][m-1] = gHat.get(i,j+1,m);
            flf_local[0][m-1]  =  fLF.get(i,j,m);
            flf_local[1][m-1]  =  fLF.get(i+1,j,m);
            flf_local[2][m-1]  =  gLF.get(i,j,m);
            flf_local[3][m-1]  =  gLF.get(i,j+1,m);
            qtmp[m-1] = q.get(i,j,m);
        }
        // plow = min(eps,plow);

        double rescale[2][2][2][2];
        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        for(int i3=0; i3<2; i3++)
        for(int i4=0; i4<2; i4++)
            rescale[i1][i2][i3][i4] = 1.0;

        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        for(int i3=0; i3<2; i3++)
        for(int i4=0; i4<2; i4++)
        {
            aa[0] = double(i1)*amin[i-1][j-1][0];
            aa[1] = double(i2)*amin[i-1][j-1][1];
            aa[2] = double(i3)*amin[i-1][j-1][2];
            aa[3] = double(i4)*amin[i-1][j-1][3];

            for( int m=1; m <= meqn; m++ )
            {
                ff[0] = aa[0]*fhat_local[0][m-1]+(1.0-aa[0])*flf_local[0][m-1];
                ff[1] = aa[1]*fhat_local[1][m-1]+(1.0-aa[1])*flf_local[1][m-1];
                ff[2] = aa[2]*fhat_local[2][m-1]+(1.0-aa[2])*flf_local[2][m-1];
                ff[3] = aa[3]*fhat_local[3][m-1]+(1.0-aa[3])*flf_local[3][m-1];
                qh[m-1] = qtmp[m-1] + dtx*(ff[0]-ff[1]) + dty*(ff[2]-ff[3]);
            }

            rho    = qh[0];
            u1     = qh[1]/rho;
            u2     = qh[2]/rho;
            u3     = qh[3]/rho;
            energy = qh[4];
            phigh  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

            if (phigh < 0.0)
            {
                root1 = 0.0;
                root2 = 1.0;

                for( int m=1; m <= meqn; m++ )
                    qtmp2[m-1] = 0.0;

                // bisection 
                // maximal iteration will affect the speed of the code
                // But reducing maximal iteration will lead to larger 
                // numerical diffusion that is introduced from low order scheme.
                for(int nn = 1; nn<5; nn++)
                {
                    rate = (root1+root2)/2.0;
                    for( int m=1; m <= meqn; m++ )
                    {
                        qtmp2[m-1] = rate*qh[m-1]+(1.0-rate)*qlf[m-1];
                    }
                    rho    = qtmp2[0];
                    u1     = qtmp2[1]/rho;
                    u2     = qtmp2[2]/rho;
                    u3     = qtmp2[3]/rho;
                    energy = qtmp2[4];
                    phigh  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));
                    if (phigh < 0.0)
                        root2 = rate;
                    else
                        root1 = rate;
                }
                if (phigh > 0.0)
                    rate = rate;
                else
                    rate = root1;

                rescale[i1][i2][i3][i4] = rate;

/*
           for( int m=1; m <= meqn; m++ )
               qtmp2[m-1] = rate*qh[m-1]+(1.0-rate)*qlf[m-1];

           rho    = qtmp2[0];
           u1     = qtmp2[1]/rho;
           u2     = qtmp2[2]/rho;
           u3     = qtmp2[3]/rho;
           energy = qtmp2[4];
           phigh2  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3));

           if (phigh2 < 0.0){
               cout << " something is wrong (i,j)="<<i<<" "<<j<<" "
                    << phigh2<<" "<< phigh<<" plow="<<plow<< endl;
               cout << aa[0]<<" "<<aa[1]<<" "
                    <<" "<<aa[2]<<" "<<aa[3]<<endl;
               cout << rate<< " "<< i1<<i2<<i3<<i4 << endl;
           }
*/

            }

        }   //end of for i1,i2 loop
        
        double rescale2[4] = {1.0, 1.0, 1.0, 1.0};
        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        for(int i3=0; i3<2; i3++)
        {
            rescale2[0] = Min(rescale2[0],rescale[1][i1][i2][i3]);
            rescale2[1] = Min(rescale2[1],rescale[i1][1][i2][i3]);
            rescale2[2] = Min(rescale2[2],rescale[i1][i2][1][i3]);
            rescale2[3] = Min(rescale2[3],rescale[i1][i2][i3][1]);
        }
            
        for(int i1=0; i1<4; i1++)
            amin[i-1][j-1][i1] = rescale2[i1]*amin[i-1][j-1][i1];
    }



    /////////////////////////////////////////////////////////
    // Share information with neighbors to determine thetas
    /////////////////////////////////////////////////////////
    for (int j=1; j<=my; j++)
    {
        for (int i=1; i<=mx-1; i++)
        {
            thex[i][j-1] = Min(amin[i][j-1][0],amin[i-1][j-1][1]);
        }
        thex[0][j-1] = amin[0][j-1][0];
        thex[mx][j-1]= amin[mx-1][j-1][1];
    }

    for (int i=1; i<=mx; i++)
    {
        for (int j=1; j<=my-1; j++)
        {
            they[i-1][j] = Min(amin[i-1][j][2],amin[i-1][j-1][3]);
        }
        they[i-1][0] = amin[i-1][0][2];
        they[i-1][my]= amin[i-1][my-1][3];
    }

    ///////////////////////////////////////////
    // Limit both fluxes (two directions each)
    ///////////////////////////////////////////
    for (int i=1; i<=mx+1; i++)
    for (int j=1; j<=my; j++)
    {
        for( int m=1; m <= meqn; m++ )
        {
            ftmp = thex[i-1][j-1]*(fHat.get(i,j,m)-fLF.get(i,j,m))+fLF.get(i,j,m);
            fHat.set(i,j,m,ftmp);
        }
    }

    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my+1; j++)
    { 
        for( int m=1; m <= meqn; m++ )
        {
            ftmp = they[i-1][j-1]*(gHat.get(i,j,m)-gLF.get(i,j,m))+gLF.get(i,j,m);
            gHat.set(i,j,m,ftmp);
        }
    }

}
