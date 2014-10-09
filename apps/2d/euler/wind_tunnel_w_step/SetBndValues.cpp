#include <cmath>
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"

void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC3& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));


// This is a user-supplied routine that sets the the boundary conditions
//
void SetBndValues( StateVars& Q )
{

    dTensorBC3& q   = Q.ref_q();
    dTensorBC3& aux = Q.ref_aux();

    const int space_order = global_ini_params.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();

    if( mx%5 != 0 || my%5 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select mx and my to be a multiple of 5.\n");
        exit(1);
    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // LEFT BOUNDARY (copy initial conditions)
    // ********************************************************************* //
    void QinitFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals);
    SampleFunction( 1-mbc, 0, 1-mbc, my+mbc, q, aux, q, &QinitFunc );
    // ********************************************************************* //

    // ********************************************************************* //
    // RIGHT BOUNDARY (Zeroth-order extrapolation)
    // ********************************************************************* //
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my; j++)
    for (int m=1; m<=meqn; m++)
    {
        double tmp = q.get(mx,j,m);                    
        q.set(i,j,m, tmp );
    }
    // ********************************************************************* //

    // Compute index where the step is located.  
    // q(istep,:) is inside the domain.
    const int istep = (mx / 5);
    const int jstep = (my / 5)+1;

    // ********************************************************************* //
    // TOP BOUNDARY
    // ********************************************************************* //
    for( int i=1-mbc; i<= mx+mbc; i++)
    for( int j=my+1; j <= (my+mbc); j++ )
    {
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(i, 2*my+1-j, m);                    
            q.set(i,j,m, tmp );

        }

        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );

    }

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 0.6 and y = 0;
    //
    //      Segment 2: x = 0.6; 0 < y < 0.2
    //
    //      Segment 3: 0.6 < x < 3.0; y = 0.2
    //
    // ********************************************************************* //

    // Segment 1
    for (int i=1; i<= istep; i++)
    for (int j=0; j>=(1-mbc); j--)
    {
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(i,1-j, m);                    
            q.set(i,j,m, tmp );

        }
        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );
    }

    // Segment 2
    for (int i = istep+1; i <= istep+mbc; i++)
    for (int j = 1; j <= jstep; j++ )
    {
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(2*istep+1-i, j, m);                    
            q.set(i,j,m, tmp );

        }
        // Flip the momemtum u1:
        q.set(i,j,2, -q.get(i,j,2) );
    }

    // Segment 3
    for( int i = istep+1; i <= mx+mbc; i++)
    for( int j = jstep-1; j >= jstep-mbc; j-- )
    {
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(i, 2*jstep-1-j, m);                    
            q.set(i,j,m, tmp );

        }
        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );
    }

}      


void SetBndValuesX( StateVars& Q )
{

    dTensorBC3& q   = Q.ref_q();
    dTensorBC3& aux = Q.ref_aux();

    const int space_order = global_ini_params.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();

    // Compute index where the step is located:
    if( mx%5 != 0 || my%5 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select mx and my to be a multiple of 5.\n");
        exit(1);
    }
    const int istep = (mx / 5);
    const int jstep = (my / 5)+1;


    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // LEFT BOUNDARY (copy initial conditions)
    // ********************************************************************* //
    void QinitFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals);
    SampleFunction( 1-mbc, 0, 1, my, q, aux, q, &QinitFunc );
    // ********************************************************************* //

    // ********************************************************************* //
    // RIGHT BOUNDARY (Zeroth-order extrapolation)
    // ********************************************************************* //
    for (int i = (mx+1); i<=(mx+mbc); i++)
    for (int j = jstep; j<=my; j++)
    for (int m=1; m<=meqn; m++)
    {
        double tmp = q.get(mx,j,m);                    
        q.set(i,j,m, tmp );
    }
    // ********************************************************************* //

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 0.6 and y = 0;
    //
    //      Segment 2: x = 0.6; 0 < y < 0.2
    //
    //      Segment 3: 0.6 < x < 3.0; y = 0.2
    //
    // ********************************************************************* //

    // Deal with the corner
    const double gamma = global_ini_params.get_gamma();
    const double   gm1 = gamma-1.0;

    const int jstepm1 = jstep-1;

    // density, u1, u2 and pressure (these get overwritten later)
    double den  = q.get(istep,jstepm1,1);
    double vex  = q.get(istep,jstepm1,2)/den;
    double vey  = q.get(istep,jstepm1,3)/den;
    double q2   = vex*vex+vey*vey;
    double pre  = gm1*( q.get(istep,jstepm1, 5) - 0.5*den*q2);

    // These two quantities are used for the entropy "fix"
    const double enth = ( q.get(istep,jstepm1,5) + pre  )/den;
    const double ent  = pre/pow( fabs(den), gamma ); 

    for (int j = 0; j <= 1; j++ )
    for (int i = istep+1; i <= istep + 4 - 2*j; i++)
    {

        // From Shu's code (check the exact index here ... )
        int jj = j + jstepm1 + 1;

        den  = q.get(i,jj,1);                       // density
        vex  = q.get(i,jj,2)/den;                   // u1
        vey  = q.get(i,jj,3)/den;                   // u2
        q2   = vex*vex+vey*vey;         
        pre  = gm1*( q.get(i,jj, 5) - 0.5*den*q2);  // pressure

        den  = pow( fabs(pre/ent), (1./gamma) );
        double qq2 = fabs(2.*(enth-gamma*pre/den/gm1));

        // Rescale velocities according to t0
        double t0 = sqrt(qq2/q2);
        vex       = vex*t0;
        vey       = vey*t0;

        q.set(i,jj,1, den     );
        q.set(i,jj,2, den*vex );
        q.set(i,jj,3, den*vey );
        q.set(i,jj,4, 0.      );
        q.set(i,jj,5, pre/gm1 + 0.5*den*qq2 );

    }


/*
 * This section is directly from Shu's code.  It is set only during the call to SetBCX.
 *
 *

c treat the singularity at the corner

      den  = u(nxmid,nymid,1)
      vex  = u(nxmid,nymid,2)/den
      vey  = u(nxmid,nymid,3)/den
      q2   = vex*vex+vey*vey
      pre  = gm1*(u(nxmid,nymid,4)-0.5*den*q2)
      enth = (u(nxmid,nymid,4)+pre)/den
      ent  = pre/abs(den)**gamma

      do 10 j=0,1
      do 10 i=nxmid+1,nxmid+4-2*j
        jj  = j+nymid+1
        den = u(i,jj,1)
        vex = u(i,jj,2)/den
        vey = u(i,jj,3)/den
        q2  = vex*vex+vey*vey
        pre = gm1*(u(i,jj,4)-0.5*den*q2)

        den=(abs(pre/ent))**(1./gamma)
        qq2=abs(2.*(enth-gamma*pre/den/gm1))
        t0=sqrt(qq2/q2) 
        vex=vex*t0
        vey=vey*t0

        u(i,jj,1)=den
        u(i,jj,2)=den*vex
        u(i,jj,3)=den*vey
        u(i,jj,4)=pre/gm1+0.5*den*qq2

10    continue
*/


    // Segment 2
    for (int i = istep+1; i <= istep+mbc; i++)
    for (int j = 1; j <= jstep; j++ )
    {
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(2*istep+1-i, j, m);                    
            q.set(i,j,m, tmp );

        }
        // Flip the momemtum u1:
        q.set(i,j,2, -q.get(i,j,2) );
    }

}

void SetBndValuesY( StateVars& Q )
{

    dTensorBC3& q   = Q.ref_q();
    dTensorBC3& aux = Q.ref_aux();

    const int space_order = global_ini_params.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();

    // Compute index where the step is located:
    if( mx%5 != 0 || my%5 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select mx and my to be a multiple of 5.\n");
        exit(1);
    }
    const int istep = (mx / 5);
    const int jstep = (my / 5)+1;
//assert_lt( (istep-0.5)*dx, 0.6 );
//assert_gt( (jstep-0.5)*dy, 0.2 );
//  printf("istep*dx = %f; istep = %d; ", istep*dx, istep  );
//  printf("jstep*dy = %f; jstep = %d\n", jstep*dy, jstep );

    // ********************************************************************* //
    // TOP BOUNDARY
    // ********************************************************************* //
    for( int i=1-mbc; i<= mx+mbc; i++)
    for( int j=my+1; j <= (my+mbc); j++ )
    {
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(i, 2*my+1-j, m);                    
            q.set(i,j,m, tmp );

        }

        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );

    }

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 0.6 and y = 0;
    //
    //      Segment 2: x = 0.6; 0 < y < 0.2
    //
    //      Segment 3: 0.6 < x < 3.0; y = 0.2
    //
    // ********************************************************************* //

    // Segment 1
    for (int i=1-mbc; i <= istep; i++)
    for (int j=0; j>=(1-mbc); j--)
    {
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(i,1-j, m);                    
            q.set(i,j,m, tmp );

        }
        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );
    }

    // Segment 3
    for( int i = istep+1; i <= mx+mbc; i++)
    for( int j = jstep-1; j >= jstep-mbc; j-- )
    {
//printf("Setting j = %d from %d\n", j, 2*jstep-1-j);
        for (int m=1; m<=meqn; m++)
        {

            double tmp = q.get(i, 2*jstep-1-j, m);                    
            q.set(i,j,m, tmp );

        }
        // Flip the momemtum u2:
        q.set(i,j,3, -q.get(i,j,3) );
    }


}
