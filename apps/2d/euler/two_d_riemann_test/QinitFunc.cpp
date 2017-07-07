#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "IniParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
//      2D Euler equations.
//
//      Q(1) = density
//      Q(2) = rho u^1
//      Q(3) = rho u^2
//      Q(4) = rho u^3
//      Q(5) = Energy
//
// Input:
//
//    xpts( 1:numpts, 1:2 )      - The x, and y-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// In this file we specify initial conditions for a number of
// 2d Riemann problems that have been studied in many papers.  See
//
//      "A posteriori subcell limiting of the discontinuous Galerkin finite
//      element method for hyperbolic conservation laws," Dumbser et al, JCP (2014).
//
// or 
//
//      "Solution of Two-Dimensional Riemann Problems for Gas Dynamics without 
//       Riemann Problem Solvers," Kurganov and Tadmor. 
//
// The Moe-Rossmanith-Seal moment limiter paper
// (https://arxiv.org/abs/1507.03024) has references and a description of RP1 
// and RP3.
//     
// We use the indicator riemann_problem_number in the [euler] section of the parameters file to 
// indicate which Riemann problem we are solving.
//
// See also: $FINESS/lib/2d/blanks/QinitFunc
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    // Gas consant
    const double gamma = global_ini_params.get_gamma();

    // Riemann problem number
    const int    riemann_number = global_ini_params.get_riemann_problem_number();

    // Loop over grid points
    const int numpts=xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {

        // Density, velocity, and pressure for each corner in the Riemann
        // problem.
        double rho1,rho2,rho3,rho4;
        double u11,u12,u13,u14;
        double u21,u22,u23,u24;
        double press1,press2,press3,press4;
 
        // Default value is to place the discontinuity at the origin.  The only
        // Riemann problem that doesn't have this is RP1 becuase it moves so
        // fast.
        //
        // (TODO - move this to a user parameter inside parameters.ini? -DS)
        double xcorner=0.0;
        double ycorner=0.0;

        // Density, velocity and pressure
        double rho,u1,u2,u3,press;
        u3 = 0.0;

        // TODO - what is this section on?  I think it makes sense to pass
        // this information in through the riemann_number paramter.
//      int center=0;       
//      if(center==0)
//      {
//          rho1=0.532258064516129;
//          u11 = 1.206045378311055;
//          u21 = 0.0;
//          u3 = 0.0;
//          press1=0.3;

//          rho2   = 0.137992831541219;
//          u12    = 1.206045378311055;
//          u22    = 1.206045378311055;
//          u3    =  0.0;
//          press2 =  0.029032258064516;

//          rho3   =  0.532258064516129;
//          u13    =  0.0;
//          u23    =  0.532258064516129;
//          u3    =  0.0;
//          press3 =  0.3;

//          rho4   = 1.5;
//          u14    = 0.0;
//          u24    = 0.0;
//          u3    = 0.0;
//          press4 = 1.5;

//          xcorner=0.3;
//          ycorner=0.3;
//      }

        if( riemann_number == 1 )
        {

            rho1=0.5323;
            u11 = 1.206;
            u21 = 0.0;
            u3 = 0.0;
            press1=0.3;

            rho2   = 0.138;
            u12    = 1.206;
            u22    = 1.206;
            u3    =  0.0;
            press2 =  0.029;

            rho3   =  0.5323;
            u13    =  0.0;
            u23    =  1.206;
            u3    =  0.0;
            press3 =  0.3;

            rho4   = 1.5;
            u14    = 0.0;
            u24    = 0.0;
            u3    = 0.0;
            press4 = 1.5;

            xcorner = 0.3;
            ycorner = 0.3;

        }

        if(riemann_number==2)
        {

            rho1=0.5065;
            u11 = 0.8939;
            u21 = 0.0;
            u3 = 0.0;
            press1=0.35;

            rho2   = 1.1;
            u12    =  0.8939;
            u22    =  0.8939;
            u3    =  0.0;
            press2 =  1.1;

            rho3   =  0.5065;
            u13    =  0.0;
            u23    =  0.8939;
            u3    =  0.0;
            press3 =  0.35;

            rho4   = 1.1;
            u14    = 0.0;
            u24    = 0.0;
            u3    = 0.0;
            press4 = 1.1;

            xcorner = 0.0;
            ycorner = 0.0;

        }
        if(riemann_number==3)
        {

            rho1    = 2.0;
            u11     = 0.75;
            u21     = 0.5;
            u3      = 0.0;
            press1  = 1.0;

            rho2    =  1.0;
            u12     = -0.75;
            u22     =  0.5;
            u3      =  0.0;
            press2  =  1.0;

            rho3    =  3.0;
            u13     = -0.75;
            u23     = -0.5;
            u3      =  0.0;
            press3  =  1.0;

            rho4   = 1.0;
            u14    = 0.75;
            u24    = -0.5;
            u3     = 0.0;
            press4 = 1.0;

            xcorner = 0.0;
            ycorner = 0.0;

        }

        if(riemann_number==4)
        {

            rho1    =  1.0;
            u11     = -0.6259;
            u21     = 0.1;
            u3      = 0.0;
            press1  = 1.0;

            rho2    = 0.8;
            u12     = 0.1;
            u22     = 0.1;
            u3      = 0.0;
            press2  = 1.0;

            rho3    = 1.0;
            u13     = 0.1;
            u23     =-0.6259;
            u3      = 0.0;
            press3  = 1.0;

            rho4    = 0.5917;
            u14     = 0.1;
            u24     = 0.1;
            u3      = 0.0;
            press4  = 0.4;

            xcorner = 0.0;
            ycorner = 0.0;

        }
        if(riemann_number==5)
        {

            rho1=1.0;
            u11 = 0.7276;
            u21 = 0.0;
            u3 = 0.0;
            press1=1.0;

            rho2   = 0.8;
            u12    =  0.0;
            u22    =  0.0;
            u3    =  0.0;
            press2 =  1.0;

            rho3   =  1.0;
            u13    =  0.0;
            u23    =  0.7276;
            u3    =  0.0;
            press3 =  1.0;

            rho4   = 0.5313;
            u14    = 0.0;
            u24    = 0.0;
            u3    = 0.0;
            press4 = 0.4;

            xcorner = 0.0;
            ycorner = 0.0;

        }

        // Determine which quadrant
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        if(x < xcorner)
        {
            if(y<ycorner)
            {
                rho   = rho2;
                u1    =  u12;
                u2    =  u22;
                u3    =  0.0;
                press =  press2;}
            else
            {
                rho=rho1;
                u1 = u11;
                u2 = u21;
                u3 = 0.0;
                press=press1;
            }
        }
        else
        {
            if(y<ycorner)
            {
                rho   =  rho3;
                u1    =  u13;
                u2    =  u23;
                u3    =  0.0;
                press =  press3;
            }
            else
            {
                rho   = rho4;
                u1    = u14;
                u2    = u24;
                u3    = 0.0;
                press = press4;
            }
             
        } 

        // Derived conserved variable
        double energy = press/(gamma-1.0e0) + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
