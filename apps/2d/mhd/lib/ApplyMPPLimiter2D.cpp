#include <cmath>
#include <iostream>
#include "stdio.h"
#include "dog_math.h"
#include "IniParams.h"
#include "constants.h"
#include "tensors.h"
using namespace std;

// MPP limiter for 2D Ideal MHD equations.
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

    // Parameters specific for MHD equations
    const double gamma = global_ini_params.get_gamma();
    const double gm1 = gamma - 1.0;

    // Storage
    dTensor2 gmin(mx, my);
    dTensor3 amin(mx, my, 4);
    dTensor2 thex(mx+1, my);
    dTensor2 they(mx, my+1);
    double ff[4];
    double aa[4];
    //double qh[meqn],qtmp[meqn],qtmp2[meqn],qtmp3[meqn],qlf[meqn],fhat_local[4][meqn],flf_local[4][meqn];
    dTensor1 qh(meqn), qtmp(meqn), qtmp2(meqn), qtmp3(meqn), qlf(meqn);
    dTensor2 fhat_local(4, meqn), flf_local(4, meqn);
    double rho,u1,u2,u3,energy, B1, B2, B3, plow,phigh,ffsum,ftmp,phigh2;
    double ac,bc,cc,delta,root1,root2,rate;
    int isgn[4];

    /////////////////////////////////////////
    // Step 1: limit rho
    /////////////////////////////////////////
    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my; j++)
        gmin.set(i, j, rho_min - ( q.get(i,j,1) 
                       + dtx*( fLF.get(i,j,1) - fLF.get(i+1,j,1))
                       + dty*( gLF.get(i,j,1) - gLF.get(i,j+1,1))));
 
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
                amin.set(i, j, k+1, Min(gmin.get(i, j)/(ffsum-eps),1.e0));
            else
                amin.set(i, j, k+1, 1.e0);
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
            qlf.set(m, q.get(i,j,m) + dtx*(fLF.get(i,j,m)-fLF.get(i+1,j,m))
                                  + dty*(gLF.get(i,j,m)-gLF.get(i,j+1,m)));

        rho    = qlf.get(1);
        u1     = qlf.get(2)/rho;
        u2     = qlf.get(3)/rho;
        u3     = qlf.get(4)/rho;
        energy = qlf.get(5);
        B1     = qlf.get(6);
        B2     = qlf.get(7);
        B3     = qlf.get(8);
        plow   = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3) - 0.5*(B1*B1 + B2*B2 + B3*B3));

        // Quick error check
        if (plow < 0.0){
	    cerr << "(" << i << ", " << j << ")  ";
	    cerr << "plow=" << plow << "  ";
          cerr << "Negative solution in Lax-Friedrich's flux" << endl;
	}

        for( int m=1; m <= meqn; m++ )
        {
            fhat_local.set(1, m, fHat.get(i,j,m));
            fhat_local.set(2, m, fHat.get(i+1,j,m));
            fhat_local.set(3, m, gHat.get(i,j,m));
            fhat_local.set(4, m, gHat.get(i,j+1,m));
            flf_local.set(1, m,  fLF.get(i,j,m));
            flf_local.set(2, m,  fLF.get(i+1,j,m));
            flf_local.set(3, m,  gLF.get(i,j,m));
            flf_local.set(4, m,  gLF.get(i,j+1,m));
            qtmp.set(m, q.get(i,j,m));
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
            aa[0] = double(i1)*amin.get(i, j, 1);
            aa[1] = double(i2)*amin.get(i, j, 2);
            aa[2] = double(i3)*amin.get(i, j, 3);
            aa[3] = double(i4)*amin.get(i, j, 4);

            for( int m=1; m <= meqn; m++ )
            {
                ff[0] = aa[0]*fhat_local.get(1, m)+(1.0-aa[0])*flf_local.get(1, m);
                ff[1] = aa[1]*fhat_local.get(2, m)+(1.0-aa[1])*flf_local.get(2, m);
                ff[2] = aa[2]*fhat_local.get(3, m)+(1.0-aa[2])*flf_local.get(3, m);
                ff[3] = aa[3]*fhat_local.get(4, m)+(1.0-aa[3])*flf_local.get(4, m);
                qh.set(m, qtmp.get(m) + dtx*(ff[0]-ff[1]) + dty*(ff[2]-ff[3]));
            }

            rho    = qh.get(1);
            u1     = qh.get(2)/rho;
            u2     = qh.get(3)/rho;
            u3     = qh.get(4)/rho;
            energy = qh.get(5);
            B1     = qh.get(6);
            B2     = qh.get(7);
            B3     = qh.get(8);
            phigh  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3) - 0.5*(B1*B1 + B2*B2 + B3*B3));

            if (phigh < 0.0)
            {
                root1 = 0.0;
                root2 = 1.0;

                for( int m=1; m <= meqn; m++ )
                    qtmp2.set(m, 0.0);

                // bisection
                for(int nn = 1; nn<10; nn++)
                {
                    rate = (root1+root2)/2.0;
                    for( int m=1; m <= meqn; m++ )
                    {
                        qtmp2.set(m,  rate*qh.get(m)+(1.0-rate)*qlf.get(m));
                    }
                    rho    = qtmp2.get(1);
                    u1     = qtmp2.get(2)/rho;
                    u2     = qtmp2.get(3)/rho;
                    u3     = qtmp2.get(4)/rho;
                    energy = qtmp2.get(5);
                    B1     = qtmp2.get(6);
                    B2     = qtmp2.get(7);
                    B3     = qtmp2.get(8);
                    phigh  = gm1*(energy-0.5*rho*(u1*u1+u2*u2+u3*u3) - 0.5*(B1*B1 + B2*B2 + B3*B3));
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
            amin.set(i, j, i1+1, rescale2[i1]*amin.get(i, j, i1+1));
    }



    /////////////////////////////////////////////////////////
    // Share information with neighbors to determine thetas
    /////////////////////////////////////////////////////////
    for (int j=1; j<=my; j++)
    {
        for (int i=1; i<=mx-1; i++)
        {
            thex.set(i+1, j, Min(amin.get(i+1, j, 1), amin.get(i, j, 2)));
        }
        thex.set(1, j, amin.get(1, j, 1));
        thex.set(mx+1, j, amin.get(mx, j, 2));
    }

    for (int i=1; i<=mx; i++)
    {
        for (int j=1; j<=my-1; j++)
        {
            they.set(i, j+1, Min(amin.get(i, j+1, 3), amin.get(i, j, 4)));
        }
        they.set(i, 1, amin.get(i, 1, 3));
        they.set(i, my+1, amin.get(i, my, 4));
    }

    ///////////////////////////////////////////
    // Limit both fluxes (two directions each)
    ///////////////////////////////////////////
    for (int i=1; i<=mx+1; i++)
    for (int j=1; j<=my; j++)
    {
        for( int m=1; m <= meqn; m++ )
        {
            ftmp = thex.get(i, j) * (fHat.get(i,j,m)-fLF.get(i,j,m))+fLF.get(i,j,m);
            fHat.set(i,j,m,ftmp);
        }
    }

    for (int i=1; i<=mx; i++)
    for (int j=1; j<=my+1; j++)
    { 
        for( int m=1; m <= meqn; m++ )
        {
            ftmp = they.get(i, j)*(gHat.get(i,j,m)-gLF.get(i,j,m))+gLF.get(i,j,m);
            gHat.set(i,j,m,ftmp);
        }
    }

}
