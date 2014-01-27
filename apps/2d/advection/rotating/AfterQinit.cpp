#include<cmath>
#include "dogdefs.h"
#include "DogParamsCart2.h"
#include "DogSolverCart2.h"

// Smooth initial conditions 
void AfterQinit(DogSolverCart2& solver)
{
    dTensorBC4& aux = solver.get_aux();
    dTensorBC4& q  = solver.get_q();
    const int   mx = q.getsize(1);
    const int   my = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(2);
    const int method1 = int((sqrt(1+8*kmax)-1)/2);  
    double minmod(double,double,double);

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    const double M=10.0;
    const double Mdx2=M*dx*dx;
    const double Mdy2=M*dy*dy;

    int m=1;
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)    
        for (int j=(2-mbc); j<=(my+mbc-1); j++)          
            for (int m=1; m<=meqn; m++)
            {
                bool limited_linear = false;

                // x-direction: limit linear term
                double wx_now = q.get(i,j,m,2);

                if(fabs(wx_now) > Mdx2) // so the peaks aren't clipped
                {
                    double dwp        =  q.get(i+1,j,m,1) - q.get(i,j,m,1);
                    double dwm        =  q.get(i,j,m,1) - q.get(i-1,j,m,1);
                    double wx_limited = minmod(wx_now, dwp/2.0, dwm/2.0);

                    if(fabs(wx_now-wx_limited)>1.0e-15)
                    {
                        limited_linear=true;
                        q.set(i,j,m,2, wx_limited );
                    }
                }

                // y-direction: limit linear term      
                double wy_now     =  q.get(i,j,m,3);
                if(fabs(wy_now) > Mdy2) // so the peaks aren't clipped
                {
                    double dwp        =  q.get(i,j+1,m,1) - q.get(i,j,m,1);
                    double dwm        =  q.get(i,j,m,1) - q.get(i,j-1,m,1);
                    double wy_limited = minmod(wy_now, dwp/2.0, dwm/2.0);

                    if(fabs(wy_now-wy_limited)>1.0e-15)
                    {
                        limited_linear = true;
                        q.set(i,j,m,3, wy_limited );
                    }
                }

                // if linear limited terms were limited,
                // then zero out all high-order terms
                if(limited_linear)
                {
                    for (int k=4; k<=kmax; k++)
                    {
                        q.set(i,j,m,k, 0.0 );
                    }
                }

            }
}
