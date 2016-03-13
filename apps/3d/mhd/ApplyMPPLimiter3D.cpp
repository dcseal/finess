#include "IniParams.h"
#include "constants.h"
#include "tensors.h"
#include "dog_math.h"

#include "iostream"

using namespace std;

// MPP limiter for 3D Ideal MHD equations.
//
// See: "http://arxiv.org/abs/1411.0328" and references therein for more
// details.
void ApplyMPPLimiter3D( 
        const double dt, const dTensorBC4& q, 
        const dTensorBC4& fLF, const dTensorBC4& gLF, const dTensorBC4& hLF,
        dTensorBC4& fHat, dTensorBC4& gHat, dTensorBC4& hHat)
{

    // Parameters for the current grid
    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int     mz = global_ini_params.get_mz();
    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();
    const double dz   = global_ini_params.get_dz();

    // Parameters used for limiting
    const double eps = 1.0e-12;
    const double dtx = dt/dx, dty = dt/dy, dtz = dt/dz;
    const double rho_min = eps;

    // Parameters specific for MHD equations
    const double gamma = global_ini_params.get_gamma();
    const double gm1 = gamma - 1.0;

    // Storage
    dTensor3 gmin(mx, my, mz);
    dTensor4 amin(mx, my, mz, 6);
    dTensor3 thex(mx+1, my, mz);
    dTensor3 they(mx, my+1, mz);
    dTensor3 thez(mx, my, mz+1);
    //double qh[meqn],qtmp[meqn],qtmp2[meqn],qtmp3[meqn],qlf[meqn],fhat_local[4][meqn],flf_local[4][meqn];

    /////////////////////////////////////////
    // Step 1: limit rho
    /////////////////////////////////////////
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
            for (int k=1; k<=mz; k++)
                gmin.set(i, j, k, rho_min - ( q.get(i,j,k,1) 
                            + dtx*( fLF.get(i,j,k,1) - fLF.get(i+1,j,k,1))
                            + dty*( gLF.get(i,j,k,1) - gLF.get(i,j+1,k,1))
                            + dtz*( hLF.get(i,j,k,1) - hLF.get(i,j,k+1,1))));
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
            for (int k=1; k<=mz; k++)
            {
                double ff[6];
                ff[0] = dtx*(fHat.get(i,j,k,1)  -fLF.get(i,j,k,1)  );
                ff[1] =-dtx*(fHat.get(i+1,j,k,1)-fLF.get(i+1,j,k,1));
                ff[2] = dty*(gHat.get(i,j,k,1)  -gLF.get(i,j,k,1)  );
                ff[3] =-dty*(gHat.get(i,j+1,k,1)-gLF.get(i,j+1,k,1));
                ff[4] = dtz*(hHat.get(i,j,k,1)  -hLF.get(i,j,k,1)  );
                ff[5] =-dtz*(hHat.get(i,j,k+1,1)-hLF.get(i,j,k+1,1));

                int isgn[6];
                for (int kk=0; kk<6; kk++)
                {
                    if (ff[kk]<0.0)
                        isgn[kk] = 1;
                    else
                        isgn[kk] = 0;
                }
                double ffsum;
                ffsum = double(isgn[0])*ff[0];

                for (int kk=1; kk<6; kk++)
                    ffsum = ffsum + double(isgn[kk])*ff[kk];

                for (int kk=0; kk<6; kk++)
                {
                    if (isgn[kk]==1)
                        amin.set(i, j, k, kk+1, Min(gmin.get(i, j, k)/(ffsum-eps),1.e0));
                    else
                        amin.set(i, j, k, kk+1, 1.e0);
                }
            }

    /////////////////////////////////////////
    // Step 2: limit pressure
    /////////////////////////////////////////
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
            for (int k=1; k<=mz; k++)
            {

                dTensor1 qh(meqn), qtmp(meqn), qtmp2(meqn),  qlf(meqn);
                // Write the update for a LF solve
                for( int m=1; m <= meqn; m++ )  
                    qlf.set(m, q.get(i,j,k,m) + dtx*(fLF.get(i,j,k,m)-fLF.get(i+1,j,k,m))
                            + dty*(gLF.get(i,j,k,m)-gLF.get(i,j+1,k,m))
                            + dtz*(hLF.get(i,j,k,m)-hLF.get(i,j,k+1,m)));
                double rho,u1,u2,u3,energy, B1, B2, B3, plow,phigh,ftmp;

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
                    cerr << "(" << i << ", " << j << ", " << k << ")  ";
                    cerr << "plow=" << plow << "  ";
                    cerr << "Negative solution in Lax-Friedrich's flux. \n";
                }
                dTensor2 fhat_local(6, meqn), flf_local(6, meqn);
                for( int m=1; m <= meqn; m++ )
                {
                    fhat_local.set(1, m, fHat.get(i,j,k,m));
                    fhat_local.set(2, m, fHat.get(i+1,j,k,m));
                    fhat_local.set(3, m, gHat.get(i,j,k,m));
                    fhat_local.set(4, m, gHat.get(i,j+1,k,m));
                    fhat_local.set(5, m, hHat.get(i,j,k,m));
                    fhat_local.set(6, m, hHat.get(i,j,k+1,m));
                    flf_local.set(1, m,  fLF.get(i,j,k,m));
                    flf_local.set(2, m,  fLF.get(i+1,j,k,m));
                    flf_local.set(3, m,  gLF.get(i,j,k,m));
                    flf_local.set(4, m,  gLF.get(i,j+1,k,m));
                    flf_local.set(5, m,  hLF.get(i,j,k,m));
                    flf_local.set(6, m,  hLF.get(i,j,k+1,m));
                    qtmp.set(m, q.get(i,j,k,m));
                }
                // plow = min(eps,plow);


                double rescale[2][2][2][2][2][2];
                for(int i1=0; i1<2; i1++)
                    for(int i2=0; i2<2; i2++)
                        for(int i3=0; i3<2; i3++)
                            for(int i4=0; i4<2; i4++)
                                for(int i5=0; i5<2; i5++)
                                    for(int i6=0; i6<2; i6++)
                                        rescale[i1][i2][i3][i4][i5][i6] = 1.0;

                double aa[6];
                for(int i1=0; i1<2; i1++)
                    for(int i2=0; i2<2; i2++)
                        for(int i3=0; i3<2; i3++)
                            for(int i4=0; i4<2; i4++)
                                for(int i5=0; i5<2; i5++)
                                    for(int i6=0; i6<2; i6++)
                                    {
                                        double ff[6];
                                        aa[0] = double(i1)*amin.get(i, j, k, 1);
                                        aa[1] = double(i2)*amin.get(i, j, k, 2);
                                        aa[2] = double(i3)*amin.get(i, j, k, 3);
                                        aa[3] = double(i4)*amin.get(i, j, k, 4);
                                        aa[4] = double(i5)*amin.get(i, j, k, 5);
                                        aa[5] = double(i6)*amin.get(i, j, k, 6);
                                        for( int m=1; m <= meqn; m++ )
                                        {
                                            ff[0] = aa[0]*fhat_local.get(1, m)+(1.0-aa[0])*flf_local.get(1, m);
                                            ff[1] = aa[1]*fhat_local.get(2, m)+(1.0-aa[1])*flf_local.get(2, m);
                                            ff[2] = aa[2]*fhat_local.get(3, m)+(1.0-aa[2])*flf_local.get(3, m);
                                            ff[3] = aa[3]*fhat_local.get(4, m)+(1.0-aa[3])*flf_local.get(4, m);
                                            ff[4] = aa[4]*fhat_local.get(5, m)+(1.0-aa[4])*flf_local.get(5, m);
                                            ff[5] = aa[5]*fhat_local.get(6, m)+(1.0-aa[5])*flf_local.get(6, m);
                                            qh.set(m, qtmp.get(m) + dtx*(ff[0]-ff[1]) + dty*(ff[2]-ff[3]) + dtz*(ff[4]-ff[5]));
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
                                        //	    cerr << "(" << i << ", " << j << ", " << k << ")  " "phigh = " << phigh << "\n";
                                        double ac,bc,cc,delta,root1,root2,rate;
                                        if (phigh < 0.0)
                                        {
                                            //		cerr << "pressure limiter triggered: " << i <<", " <<j <<", " <<k << "\n";
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

                                            rescale[i1][i2][i3][i4][i5][i6] = rate;

                                        }

                                    }        

                double rescale2[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
                for(int i1=0; i1<2; i1++)
                    for(int i2=0; i2<2; i2++)
                        for(int i3=0; i3<2; i3++)
                            for(int i4=0; i4<2; i4++)
                                for(int i5=0; i5<2; i5++)
                                {
                                    rescale2[0] = Min(rescale2[0],rescale[1][i1][i2][i3][i4][i5]);
                                    rescale2[1] = Min(rescale2[1],rescale[i1][1][i2][i3][i4][i5]);
                                    rescale2[2] = Min(rescale2[2],rescale[i1][i2][1][i3][i4][i5]);
                                    rescale2[3] = Min(rescale2[3],rescale[i1][i2][i3][1][i4][i5]);
                                    rescale2[4] = Min(rescale2[4],rescale[i1][i2][i3][i4][1][i5]);
                                    rescale2[5] = Min(rescale2[5],rescale[i1][i2][i3][i4][i5][1]);
                                }

                for(int i1=0; i1<6; i1++)
                    amin.set(i, j, k, i1+1, rescale2[i1]*amin.get(i, j, k, i1+1));
            }



    /////////////////////////////////////////////////////////
    // Share information with neighbors to determine thetas
    /////////////////////////////////////////////////////////
#pragma omp parallel for
    for (int j=1; j<=my; j++)
        for (int k=1; k<=mz; k++)
        {
            for (int i=1; i<=mx-1; i++)
            {
                thex.set(i+1, j, k, Min(amin.get(i+1, j, k, 1), amin.get(i, j, k, 2)));
            }
            thex.set(1, j, k, amin.get(1, j, k, 1));
            thex.set(mx+1, j, k, amin.get(mx, j, k, 2));
        }
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int k=1; k<=mz; k++)	
        {
            for (int j=1; j<=my-1; j++)
            {
                they.set(i, j+1, k, Min(amin.get(i, j+1, k, 3), amin.get(i, j, k, 4)));
            }
            they.set(i, 1, k, amin.get(i, 1, k, 3));
            they.set(i, my+1, k, amin.get(i, my, k, 4));
        }
#pragma omp parallel for
    for(int i = 1; i <= mx; i++)
        for(int j = 1; j <= my; j++){
            for(int k = 1; k <= mz - 1; k++)
                thez.set(i, j, k+1, Min(amin.get(i, j, k+1, 5), amin.get(i, j, k, 6)));
            thez.set(i, j, 1, amin.get(i, j, 1, 5));
            thez.set(i, j, mz+1, amin.get(i, j, mz, 6));
        }
#pragma omp parallel for
    for (int i=1; i<=mx+1; i++)
        for (int j=1; j<=my; j++)
            for (int k=1; k<=mz; k++)
            {
                for( int m=1; m <= meqn; m++ )
                {
                    double the = thex.get(i, j, k);
                    double ftmp = the* (fHat.get(i,j,k,m)-fLF.get(i,j,k,m))+fLF.get(i,j,k,m);
                    fHat.set(i,j,k,m,ftmp);
                }
            }
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int k=1; k<=mz; k++)
            for (int j=1; j<=my+1; j++)
            { 
                for( int m=1; m <= meqn; m++ )
                {
                    double the = they.get(i, j, k);
                    double ftmp = the*(gHat.get(i,j,k,m)-gLF.get(i,j,k,m))+gLF.get(i,j,k,m);
                    gHat.set(i,j,k,m,ftmp);
                }
            }
#pragma omp parallel for
    for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
            for (int k=1; k<=mz+1; k++)
            {
                for( int m=1; m <= meqn; m++ )
                {
                    double the = thez.get(i, j, k);
                    double ftmp = the * (hHat.get(i,j,k,m)-hLF.get(i,j,k,m))+hLF.get(i,j,k,m);
                    hHat.set(i,j,k,m,ftmp);
                }
            }


}
