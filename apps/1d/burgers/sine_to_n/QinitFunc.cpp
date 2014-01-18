#include <cmath>
#include "constants.h"
#include "tensors.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Note: for Burger's equation with these initial conditions, the solution
// stays smooth until 1 / (2*pi) \approx 1.591549430918953e-01.
//
// See: LeVeque pg 224, "Finite Volume Methods for Hyperbolic Problems."
//
// The exact 'breaking time' is given by, 
//
//      T_b = -1 / min_x ( f''(q(x,0) ) q_x( x, 0 ) ).
//
void QinitFunc(const dTensor1& xpts, 
        dTensor2& qvals)
{
    const int numpts=xpts.getsize();

    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        qvals.set(i,1, 0.5 + sin(2.0*pi*x) );
    }
}
