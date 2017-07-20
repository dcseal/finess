#include <cmath>
#include "dogdefs.h"
#include "IniParams.h"
#include "StateVars.h"

// This is a user-supplied routine that sets the the boundary conditions
//
// The double mach reflection problem has a time-dependent boundary condition.
// The top part of the domain tracks the front of the shock, whereas the
// bottom part of the domain is supposed to be a reflective boundary condition
// along the wedge, and saved as the initial conditions along the cells that
// are in front of the wedge.
// 
void SetBndValues( StateVars& Q )
{

    dTensorBC3& q    = Q.ref_q();
    dTensorBC3& aux  = Q.ref_aux();
    const double t   = Q.get_t();

    void SampleFunction( 
        int istart, int iend,
        int jstart, int jend,
        const dTensorBC3& qin, 
        const dTensorBC3& auxin,  
              dTensorBC3& Fout,
        void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));

    const int space_order = global_ini_params.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);  

    // Grid parameters for this problem
    const double dx   = global_ini_params.get_dx();
    const double xlow = global_ini_params.get_xlow();
    const double x0   = global_ini_params.get_x0();

    // Gas constant
    const double gamma = global_ini_params.get_gamma();

    // ----------------------- //
    // BOUNDARY CONDITION LOOP
    // ----------------------- //

    // ********************************************************************* //
    // LEFT BOUNDARY (copy initial conditions)
    // ********************************************************************* //
    void LeftFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals);
    SampleFunction( 1-mbc, 0, 1-mbc, my+mbc, q, aux, q, &LeftFunc );
    // ********************************************************************* //


    // ********************************************************************* //
    // RIGHT BOUNDARY (Zeroth-order extrapolation)
    // ********************************************************************* //
    for (int i=(mx+1); i<=(mx+mbc); i++)
    for (int j=1; j<=my+mbc; j++)
    for (int m=1; m<=meqn; m++)
    {
        double tmp = q.get(mx,j,m);                    
        q.set(i,j,m, tmp );
    }
    // ********************************************************************* //


    // ********************************************************************* //
    // BOTTOM BOUNDARY (Identical to initial condtitions to the left of x0)
    // ********************************************************************* //
    void BotFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
            dTensor2& qvals);
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, 0, q, aux, q, &BotFunc );

    // Make q an even function to the right of x0 along the bottom boundary
    for (int i=1; i<=mx; i++)
    for (int j=0; j>=(1-mbc); j--)
    {

        double x = xlow + (double(i)-0.5)*dx;      
        if( x >= (x0-1.0e-12) )
        {
            for (int m=1; m<=meqn; m++)
            {

                double tmp = q.get(i,1-j, m);                    
                q.set(i,j,m, tmp );

            }

            // TODO - NEED TO HAVE DIFFERENT BOUNDARY CONDITIONS APPLIED TO F, G THAN Q!
            double tmp = -q.get(i,j,3);
            q.set(i,j,3, tmp );

        }

    }
    // ********************************************************************* //

    // ***********************************************
    // TOP BOUNDARY (This part is time dependent!)
    //************************************************  
//  void TopFunc(const dTensor2& xpts,
//      const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
//          dTensor2& qvals);
//  SampleFunction( 1-mbc, mx+mbc, my+1, my+mbc, q, aux, q, &TopFunc );
    for(int i=1-mbc; i<=mx+mbc; i++)
    for(int j=my+1;  j<=my+mbc; j++)
    {

        double x = xlow + (double(i)-0.5)*dx;      
        double rho,u1,u2,u3,press;
        if ( x < x0+(20.0*t+1.0)*osq3)
        {
            rho   =  8.0;
            u1    =  8.25*cos(pi/6.0);
            u2    = -8.25*sin(pi/6.0);
            u3    =  0.0;
            press =  116.5;
        }
        else
        {
            rho   =  1.4;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            press =  1.0;
        }	

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        q.set(i,j,1, rho    );
        q.set(i,j,2, rho*u1 );
        q.set(i,j,3, rho*u2 );
        q.set(i,j,4, rho*u3 );
        q.set(i,j,5, energy );      
    }

// printf("t = %2.3e\n", t );

}




// This function is idential to the "left hand" initial conditions of the
// problem.
void LeftFunc(const dTensor2& xpts,
    const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals)
{
    const int numpts   = xpts.getsize(1);
    const double gamma = global_ini_params.get_gamma();

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double rho,u1,u2,u3,press;

        rho   =  8.0;
        u1    =  8.25*cos(pi/6.0);
        u2    = -8.25*sin(pi/6.0);
        u3    =  0.0;
        press =  116.5;

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );      
    }
}

// We use the same "initial conditions"
void BotFunc(const dTensor2& xpts,
    const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals)
{
    const int numpts   = xpts.getsize(1);
    const double gamma = global_ini_params.get_gamma();
    const double x0    = global_ini_params.get_x0();

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double rho,u1,u2,u3,press;

        if ( x < x0 )
        {
            rho   =  8.0;
            u1    =  8.25*cos(pi/6.0);
            u2    = -8.25*sin(pi/6.0);
            u3    =  0.0;
            press =  116.5;
        }
        else
        {
            rho   =  1.4;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            press =  1.0;
        }	

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );      
    }
}

// Wrappers for main Euler library
void SetBndValuesX(StateVars& Q)
{ SetBndValues( Q ); }

void SetBndValuesY(StateVars& Q)
{ SetBndValues( Q ); }


