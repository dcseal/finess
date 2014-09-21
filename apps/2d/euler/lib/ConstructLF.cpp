#include <iostream>
#include "dogdefs.h"
#include "dog_math.h"
#include "stdlib.h"
#include "IniParams.h"       


using std::cout;
using std::endl;
// User supplied functions defining the Flux function, Jacobian, and
// Hessian of the flux function.
void FluxFunc(const dTensor2&,const dTensor2&,const dTensor2&,dTensor3&);

// Used for construcing the flux function
void SampleFunction( 
    int istart, int iend,
    int jstart, int jend,
    const dTensorBC3& qin, 
    const dTensorBC3& auxin,  
          dTensorBC4& Fout,
    void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor3&));

// First-order Lax-Friedrichs solver where the boundary conditions are applied
// twice - once for constructing F, and second for constructing G.  This is used
// for the problems that have "geometry".
//
// See also ConstructLF. (same file)
void ConstructLF2( const double dt, 
    dTensorBC3& aux, dTensorBC3& q,
    dTensorBC3& F, dTensorBC3& G,
    dTensorBC3& Lstar,
    dTensorBC3& smax)

{    
    void SetBndValuesX(dTensorBC3& aux, dTensorBC3& q); // Only set conditions along x-direction
    void SetBndValuesY(dTensorBC3& aux, dTensorBC3& q); // Only set conditions along y-direction

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();

    dTensorBC3  Fhat(mx+1, my,   meqn, mbc );
    dTensorBC3  Ghat(mx,   my+1, meqn, mbc );

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC4 R( mx, my, meqn, 2, mbc );  // place-holder for the flux function

const int mbc_small      = 3;

    // Compute finite difference approximations on all of the conserved
    // variables:
    SetBndValuesX(aux, q);
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, R, &FluxFunc );

    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    for( int j = 1-mbc_small; j <= my+mbc_small; j++ )
    {
        // FINAL STEP: save the time-integrated values:
        for( int m=1; m<=meqn; m++ )
        {
            F.set(i,j,m, R.get(i,j,m,1) );
        }
    }

    double alpha_x = 1.0e-15;
    const double gamma = global_ini_params.get_gamma();

    for (int i = 1; i <= mx;  i++)
    for (int j = 1; j <= my;  j++)
    {
        // -- Compute a local wave speed -x- //
        dTensor1 xedge(2);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        const double rho = q.get(i,j,1);
        const double u1 =  q.get(i,j,2)/rho;
        const double u2 =  q.get(i,j,3)/rho;
        const double u3 =  q.get(i,j,4)/rho;
        const double energy = q.get(i,j,5);

        const double press = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

        if (press < 0.0)
            cout << " x negative press = " << press << " (x,y)= "<< 
            xedge.get(1)<<" "<<xedge.get(2) << endl;

        const double c = sqrt(fabs(gamma*press/rho));
        const double alpha = abs(u1) + c;

        smax.set( i, j, 1, Max( smax.get(i,j,1), alpha )  );

        alpha_x = Max(alpha_x, alpha);
    }

    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
        for( int m=1; m <= meqn; m++ )
        {
            double hf =  0.5*(F.get(i-1,j,m)+F.get(i,j,m)) 
                        + 0.5*alpha_x*(q.get(i-1,j,m)-q.get(i,j,m));
            Fhat.set(i,j,m, hf );
        }


    SetBndValuesY(aux, q);
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, R, &FluxFunc );

    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    for( int j = 1-mbc_small; j <= my+mbc_small; j++ )
    {
        // FINAL STEP: save the time-integrated values:
        for( int m=1; m<=meqn; m++ )
        {
            G.set(i,j,m, R.get(i,j,m,2) );        
        }
    }

    double alpha_y = 1.0e-15;


    for (int i = 1; i <= mx;  i++)
    for (int j = 1; j <= my;  j++)
    {
        // -- Compute a local wave speed -y- //
        dTensor1 xedge(2); 
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        const double rho = q.get(i,j,1);
        const double u1 =  q.get(i,j,2)/rho;
        const double u2 =  q.get(i,j,3)/rho;
        const double u3 =  q.get(i,j,4)/rho;
        const double energy = q.get(i,j,5);

        const double press = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

        if (press < 0.0)
            cout << " y negative press = " << press << " (x,y)= "<< 
            xedge.get(1)<<" "<<xedge.get(2) << endl;

  // Sound speeds
        const double c = sqrt(fabs(gamma*press/rho));
        const double alpha = abs(u2) + c;

        smax.set( i, j, 2, Max( smax.get(i,j,2), alpha )  );
        alpha_y = Max(alpha_y, alpha);
    }

    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
        for( int m=1; m <= meqn; m++ )
        {
            double hf =  0.5*(G.get(i,j-1,m)+G.get(i,j,m)) 
                        + 0.5*alpha_y*(q.get(i,j-1,m)-q.get(i,j,m));
            Ghat.set(i,j,m, hf );
        }


    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_{i,j} = Lstar = -1/dx( fh_{i+1/2,j} - fh_{i-1/2,j} )
    //                           -1/dy( gh_{i,j+1/2} - gh_{i,j-1/2} ).
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {
        printf("Error: source-term not implemented for Lax-Wendroff method\n");
        exit(1);
    }
    else  // No source term
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1,j,  m) - Fhat.get(i,j,m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,  j+1,m) - Ghat.get(i,j,m) ) / dy;
                Lstar.set(i,j,m, tmp );
            }
        }
    }

}


// First-order Lax-Friedrich's solver.  Boundary condtiions are applied once, at
// the beginning of the call to this routine.
//
// See also ConstructLF2 (same file).
void ConstructLF( const double dt, 
    dTensorBC3& aux, dTensorBC3& q,
    dTensorBC3& smax, 
    dTensorBC3& F, dTensorBC3& G,
    dTensorBC3& Lstar )
{
    void SetBndValues(dTensorBC3& aux, dTensorBC3& q); 

    // Grid and problem information
    const int mx     = global_ini_params.get_mx();
    const int my     = global_ini_params.get_my();
    const int meqn   = global_ini_params.get_meqn();
    const int maux   = global_ini_params.get_maux();
    const int mbc    = global_ini_params.get_mbc();

    const double     dx = global_ini_params.get_dx();
    const double     dy = global_ini_params.get_dy();
    const double   xlow = global_ini_params.get_xlow();
    const double   ylow = global_ini_params.get_ylow();

    dTensorBC3  Fhat(mx+1, my,   meqn, mbc );
    dTensorBC3  Ghat(mx,   my+1, meqn, mbc );

    // Sample the flux function on the entire domain:
    //
    // If "1st-order" (Euler step), then this completes this function call.
    //
    dTensorBC4 R( mx, my, meqn, 2, mbc );  // place-holder for the flux function

const int mbc_small      = 3;

    // Compute finite difference approximations on all of the conserved
    // variables:
    SetBndValues(aux, q);
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, q, aux, R, &FluxFunc );

    for( int i = 1-mbc_small; i <= mx+mbc_small; i++ )
    for( int j = 1-mbc_small; j <= my+mbc_small; j++ )
    {
        // FINAL STEP: save the time-integrated values:
        for( int m=1; m<=meqn; m++ )
        {
            F.set(i,j,m, R.get(i,j,m,1) );
            G.set(i,j,m, R.get(i,j,m,2) ); 
        }
    }

    double alpha_x = 1.0e-15, alpha_y = 1.0e-15;
    const double gamma = global_ini_params.get_gamma();

    for (int i = 1; i <= mx;  i++)
    for (int j = 1; j <= my;  j++)
    {
        // -- Compute a local wave speed //
        dTensor1 xedge(2);
        xedge.set( 1, xlow + double(i)*dx - 0.5*dx );
        xedge.set( 2, ylow + double(j)*dy - 0.5*dy );

        const double rho = q.get(i,j,1);
        const double u1 =  q.get(i,j,2)/rho;
        const double u2 =  q.get(i,j,3)/rho;
        const double u3 =  q.get(i,j,4)/rho;
        const double energy = q.get(i,j,5);

        const double press = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

        if (press < 0.0)
            cout << " x negative press = " << press << " (x,y)= "<< 
            xedge.get(1)<<" "<<xedge.get(2) << endl;

        const double c = sqrt(fabs(gamma*press/rho));
        const double alpha_xl = abs(u1) + c;
        const double alpha_yl = abs(u2) + c;

        smax.set( i, j, 1, Max( smax.get(i,j,1), alpha_xl )  );
        smax.set( i, j, 2, Max( smax.get(i,j,2), alpha_yl )  );

        alpha_x = Max(alpha_x, alpha_xl);
        alpha_y = Max(alpha_y, alpha_yl);
    }


    for (int i = 1; i <= mx+1; i++)
    for (int j = 1; j <= my;   j++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =  0.5*(F.get(i-1,j,m)+F.get(i,j,m)) 
                    + 0.5*alpha_x*(q.get(i-1,j,m)-q.get(i,j,m));
        Fhat.set(i,j,m, hf );
    }

    for (int i = 1; i<= mx;   i++)
    for (int j = 1; j<= my+1; j++)
    for( int m=1; m <= meqn; m++ )
    {
        double hf =  0.5*(G.get(i,j-1,m)+G.get(i,j,m)) 
                    + 0.5*alpha_y*(q.get(i,j-1,m)-q.get(i,j,m));
        Ghat.set(i,j,m, hf );
    }


    // --------------------------------------------------------------------- //
    // Construct Lstar, defined by:
    //
    //    d/dt q_{i,j} = Lstar = -1/dx( fh_{i+1/2,j} - fh_{i-1/2,j} )
    //                           -1/dy( gh_{i,j+1/2} - gh_{i,j-1/2} ).
    //
    // TODO - We should be able to avoid this for loop if we save Lstar in the
    // above loop without executing a second loop.  However, this requires 
    // larger strides.  (-DS)
    // --------------------------------------------------------------------- //
    if( global_ini_params.get_source_term() )
    {
        printf("Error: source-term not implemented for Lax-Wendroff method\n");
        exit(1);
    }
    else  // No source term
    {
#pragma omp parallel for
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {
            for (int m=1; m<=meqn; m++)
            {
                double tmp = -(Fhat.get(i+1,j,  m) - Fhat.get(i,j,m) ) / dx;
                tmp =  tmp   -(Ghat.get(i,  j+1,m) - Ghat.get(i,j,m) ) / dy;
                Lstar.set(i,j,m, tmp );
            }
        }
    }
}
