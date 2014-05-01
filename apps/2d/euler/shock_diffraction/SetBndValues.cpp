#include <cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "EulerParams.h"  

// This is a user-supplied routine that sets the the boundary conditions
//
void SetBndValues( dTensorBC3& aux, dTensorBC3& q )
{

    void SampleFunction( 
        int istart, int iend,
        int jstart, int jend,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  
              dTensorBC3& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

    const int space_order = dogParams.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx   = dogParamsCart2.get_dx();
    const double dy   = dogParamsCart2.get_dy();

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // RIGHT BOUNDARY (copy initial conditions)
    // ********************************************************************* //
    void QinitFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals);
    SampleFunction( mx+1, mx+mbc, 1-mbc, my+mbc, q, aux, q, &QinitFunc );
    // ********************************************************************* //

    // ********************************************************************* //
    // TOP BOUNDARY "outflow" - zeroth-order extrapolation
    // ********************************************************************* //
    for (int i=1-mbc; i<= mx+mbc; i++)
    for (int j=my+1; j<=my+mbc; j++)
    for (int m=1; m<=meqn; m++)
    {
        double tmp = q.get(i,my,m);                    
        q.set(i,j,m, tmp );
    }
    // ********************************************************************* //

    // Compute index where the step is located.  
    // q(istep,:) is inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

//assert_gt( (istep-0.5)*dx, 1.0 );
//assert_gt( (jstep-0.5)*dy, 6.0 );
//printf("xlow + (istep-0.5)*dx = %f; istep = %d; ", dogParamsCart2.get_xlow() + (istep-0.5)*dx, istep  );
//printf("ylow + (jstep-0.5)*dy = %f; jstep = %d\n", dogParamsCart2.get_ylow() + (jstep-0.5)*dy, jstep );

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // ********************************************************************* //

    // Segment 1
    for( int i = 1;       i <= istep-1; i++)
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

    // Segment 2
    for (int i = istep-1; i >= istep-mbc; i--)
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
    for (int i=istep; i<= mx+mbc; i++)
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

    // ********************************************************************* //
    // LEFT BOUNDARY (the part above the wall)
    //
    //     We consider the section x=0, 6 < y < 11 here.
    //     the other part of the "left" boundary has already been dealt with.
    // ********************************************************************* //
    SampleFunction( 1-mbc, 0, jstep, my+mbc, q, aux, q, &QinitFunc );

}      

void SetBndValuesX( dTensorBC3& aux, dTensorBC3& q )
{

    void SampleFunction( 
        int istart, int iend,
        int jstart, int jend,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  
              dTensorBC3& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

    const int space_order = dogParams.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx   = dogParamsCart2.get_dx();
    const double dy   = dogParamsCart2.get_dy();

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // RIGHT BOUNDARY (copy initial conditions)
    // ********************************************************************* //
    void QinitFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals);
    SampleFunction( mx+1, mx+mbc, 1-mbc, my+mbc, q, aux, q, &QinitFunc );
    // ********************************************************************* //

    // Compute index where the step is located.  
    // q(istep,:) is inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // ********************************************************************* //

    // Segment 2
    for (int i = istep-1; i >= istep-mbc; i--)
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

    // ********************************************************************* //
    // LEFT BOUNDARY (the part above the wall)
    //
    //     We consider the section x=0, 6 < y < 11 here.
    //     the other part of the "left" boundary has already been dealt with.
    // ********************************************************************* //
    SampleFunction( 1-mbc, 0, jstep, my+mbc, q, aux, q, &QinitFunc );

}

void SetBndValuesY( dTensorBC3& aux, dTensorBC3& q )
{

    void SampleFunction( 
        int istart, int iend,
        int jstart, int jend,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  
              dTensorBC3& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

    const int space_order = dogParams.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    const double dx   = dogParamsCart2.get_dx();
    const double dy   = dogParamsCart2.get_dy();

    if( mx%13 != 0 || my%11 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select ");
        printf("mx to be a multiple of 13, and my to be a multiple of 11.\n");
        exit(1);
    }

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // TOP BOUNDARY "outflow" - zeroth-order extrapolation
    // ********************************************************************* //
    for (int i=1-mbc; i<= mx+mbc; i++)
    for (int j=my+1; j<=my+mbc; j++)
    for (int m=1; m<=meqn; m++)
    {
        double tmp = q.get(i,my,m);                    
        q.set(i,j,m, tmp );
    }
    // ********************************************************************* //

    // Compute index where the step is located.  
    // q(istep,:) is inside the domain.
    const int istep = (mx / 13)+1;
    const int jstep = (6*(my / 11))+1;

    // ********************************************************************* //
    // BOTTOM BOUNDARY
    //
    // This region is broken into a total of three segments:
    //
    //      Segment 1: 0 < x < 1.0 and y = 6;
    //
    //      Segment 2: x = 1.0; 0 < y < 6.0;
    //
    //      Segment 3: 1.0 < x < 13.0; y = 0.0
    //
    // ********************************************************************* //

    // Segment 1
    for( int i = 1;       i <= istep-1; i++)
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

    // Segment 3
    for (int i=istep; i<= mx+mbc; i++)
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

}
