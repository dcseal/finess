#include <cmath>
#include <iostream>
#include "dogdefs.h"
#include "dog_math.h"
#include "IniParams.h"
#include "tensors.h"
#include "stdio.h"
using namespace std;

// MPP limiter for 1D Euler equations.
//
// This limiter works by considering a high order flux, fHat, and
// low-order flux, fLF, and then creates a new flux of the form:
//
//     fhat = theta fhat + (1-theta) fLF
//
// so that the conservative update, 
//
// q = q - dt/dx( f_{i+1/2} - f_{i-1/2} )
//
// retains positivity of the density and pressure.
//
// See: "http://arxiv.org/abs/1411.0328" and references therein for more
// details.
void ApplyMPPLimiter1D( const double dt, const dTensorBC2& q, const dTensorBC2& fLF, dTensorBC2& fHat )
{

    // Parameters for the current grid
    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);

    const double dx  = global_ini_params.get_dx();
    const double eps = 1.0e-12;     // TODO - do we want this hard-coded? -DS
    const double dtx = dt/dx;
    const double rho_min = eps;

    const double gamma = global_ini_params.get_gamma();
    const double gm1 = gamma - 1.0;

    // Storage
    double gmin[mx],amin[mx][2],thex[mx+1],ff[2],aa[2];
    double qh[meqn],qtmp[meqn],qlf[meqn],fhat_local[2][meqn],flf_local[2][meqn];
    double rho,u1,u2,u3,energy,plow,phigh,ffsum,ftmp;
    double ac,bc,cc,delta,root1,root2,rate;
    int isgn[2];

    // limiting on rho
    for (int i=1; i<=mx; i++)
    { gmin[i-1] = rho_min-( q.get(i,1) + dtx*( fLF.get(i,1) - fLF.get(i+1,1))); }

    for (int i=1; i<=mx; i++)
    {
        ff[0] = dtx*(fHat.get(i,1)  -fLF.get(i,1)  );
        ff[1] =-dtx*(fHat.get(i+1,1)-fLF.get(i+1,1));

        for (int k=0; k<=1; k++)
        {
            if (ff[k]<0.0)
                isgn[k] = 1;
            else
                isgn[k] = 0;
        }

        ffsum = isgn[0]*ff[0]+isgn[1]*ff[1];

        for (int k=0; k<=1; k++)
        {
            if (isgn[k]==1)
                amin[i-1][k] = Min(gmin[i-1]/(ffsum-eps),1.0);
            else
                amin[i-1][k] = 1.0;
        }
    }

    // limiting on pressure
    for (int i=1; i<=mx; i++)
    {
        for( int m=1; m <= meqn; m++ )  
            qlf[m-1]=q.get(i,m) + dtx*(fLF.get(i,m)-fLF.get(i+1,m));

        rho    = qlf[0];
        u1     = qlf[1]/rho;
        u2     = 0.0;
        u3     = 0.0;
        plow   = rho; 
        energy = plow/gm1+0.5*rho*(u1*u1);

        if (plow < 0.0)
        { cout << "Negative solution in Lax-Friedrich's flux" << endl; }

        for( int m=1; m <= meqn; m++ )
        {
            fhat_local[0][m-1] = fHat.get(i,m);
            fhat_local[1][m-1] = fHat.get(i+1,m);
            flf_local[0][m-1]  =  fLF.get(i,m);
            flf_local[1][m-1]  =  fLF.get(i+1,m);
            qtmp[m-1] = q.get(i,m);
        }
        plow = min(eps,plow);

        double rescale[2][2] = {{1.0,1.0},{1.0,1.0}};
        for(int i1=0; i1<2; i1++)
        for(int i2=0; i2<2; i2++)
        {

            aa[0] = i1*amin[i-1][0];
            aa[1] = i2*amin[i-1][1];
            for( int m=1; m <= meqn; m++ )
            {
                ff[0] = aa[0]*fhat_local[0][m-1]+(1.0-aa[0])*flf_local[0][m-1];
                ff[1] = aa[1]*fhat_local[1][m-1]+(1.0-aa[1])*flf_local[1][m-1];
                qh[m-1] = qtmp[m-1] + dtx*(ff[0]-ff[1]);
            }

            rho    = qh[0];
            u1     = qh[1]/rho;
            u2     = 0.0;
            u3     = 0.0;
            phigh  = rho; 
            energy = phigh/gm1+0.5*rho*u1*u1;

            if (phigh < 0.0)
            {
                ac = (qh[1]-qlf[1])*(qh[5]-qlf[5])-0.5*( pow(qh[2]-qlf[2],2)+pow(qh[3]-qlf[3],2)+pow(qh[4]-qlf[4],2) );
                bc = (qh[5]-qlf[5])*qlf[1]+(qh[1]-qlf[1])*qlf[5]-(qh[2]-qlf[2])*qlf[2]-(qh[3]-qlf[3])*qlf[3]-(qh[4]-qlf[4])*qlf[4]-plow/gm1*(qh[1]-qlf[1]);
                cc = qlf[1]*qlf[5]-0.5*(pow(qlf[2],2)+pow(qlf[3],2)+pow(qlf[4],2)) - plow/gm1*qlf[1];

                if (fabs(ac) >= eps)
                {
                    delta = sqrt(fabs(bc*bc-4.0*ac*cc));
                    root1 = (-bc-delta)/(2.0*ac);
                    root2 = (-bc+delta)/(2.0*ac);
                    if (ac > 0.0) 
                    {
                        rate = 1.0;
                        if (root1 >=0.0)
                            rate = Min(root1,rate);
                        if (root2 >=0.0)
                            rate = Min(root2,rate); 
                    }
                    else 
                    {
                        rate = 0.0;
                        if (root1 <=1.0)
                            rate = Max(rate,root1);
                        if (root2 <=1.0)
                            rate = Max(rate,root2); 
                    }
                }
                else 
                {
                    rate = 0.0;
                    if ( fabs(bc)>0.0 )
                    {
                        root1 = -cc/bc; 
                        if ((root1>=0.0) && (root1<=1.0))
                            rate = root1;
                    }
                }
                rescale[i1][i2] = rate;
            }

        }   //end of for i1,i2 loop
        
        double rescale2[2] = {1.0, 1.0};
        for(int i1=0; i1<2; i1++)
        {
            rescale2[0] = Min(rescale2[0],rescale[1][i1]);
            rescale2[1] = Min(rescale2[1],rescale[i1][1]);
        }
            
        for(int i1=0; i1<2; i1++)
            amin[i-1][i1] = rescale2[i1]*amin[i-1][i1];    
    }


    for (int i=1; i<=mx-1; i++)
    {
        thex[i] = Min( amin[i][0], amin[i-1][1] );
    }
    thex[0]  = amin[0][0];
    thex[mx] = amin[mx-1][1];

    // Overwrite current flux with the MPP limiter
    for (int i=1; i<=mx+1; i++)
    for( int m=1; m <= meqn; m++ )
    {

// Check whether or not we actually limited the fluxes.
//  if( fabs( 1.0-thex[i-1] ) > 1e-12 )
//  {
//  printf("thex[%d] = %2.13e\n", i-1, thex[i-1] );
//  thex[i-1] = 0.0;
//  }

        ftmp = thex[i-1]*(fHat.get(i,m)-fLF.get(i,m))+fLF.get(i,m);
        fHat.set(i,m,ftmp);
    }

}
