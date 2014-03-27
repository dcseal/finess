#include <cmath>
#include "dogdefs.h"

// This is a user-required routine.
//
// This routine defines the initial conditions for the problem.
//
// Input:
//
//    xpts( 1:numpts, 1:ndim )
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )
//
// See also: ...
//
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize();
    const double width = 2.0*0.25;

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);

        // Original function
        qvals.set( i, 1, exp( cos(2.0*pi*x) ) );

        // Derivatives of original function:
        qvals.set( i, 2, -2.0*pi*exp(cos(2.0*pi*x))*sin(2.*pi*x) );
        qvals.set( i, 3, 4.*pi*pi*(-pow(cos(2.*pi*x),2) - cos(2.*pi*x) + 1.)*exp(cos(2.*pi*x))  );

        if( fabs( x-0.5 ) < 0.2 )
        {
            qvals.set(i,3, 1.0 );
        }
        else
        {
            qvals.set(i,3,0.0);
        }
    }
  
}
