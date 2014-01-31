#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
//     ** Euler Equations **
//
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
// Inputs:
//
//     xpts( 1:numpts )         - a list of x-values at various spatial points.
//        Q( 1:numpts, 1:meqn ) - a vector of conserved variables
//      Aux( 1:numpts, 1:maux ) - vector of auxilary values
//   
// Output:
//
//    D2flux( 1:numpts, 1:meqn, 1:meqn, 1:meqn ) - f''(q) at each point.
//
// See also: FluxFunc and DFluxFunc.
void D2FluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux,
        dTensor4& D2flux)
{

    // Gas constant
    const double gamma = eulerParams.gamma;
    const double gm1   = gamma - 1.0;


    D2flux.setall(0.);  // most terms are zero

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);      

        // (Characteristic) Variables
        const double rho    = Q.get(i,1);
        const double u1     = Q.get(i,2)/rho;
        const double u2     = Q.get(i,3)/rho;
        const double u3     = Q.get(i,4)/rho;
        const double energy = Q.get(i,5);
        const double press  = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

        const double umag2 = u1*u1 + u2*u2 + u3*u3;

        // Derivatives with respect to f2:
        D2flux.set( i, 2, 1, 1, -(gm1*umag2 - 2.0*u2*u2) / rho );
        D2flux.set( i, 2, 2, 1, u1*( gamma-3.0 ) / rho );        D2flux.set( i, 2, 1, 2, u1*( gamma-3.0 ) / rho );
        D2flux.set( i, 2, 2, 2, (3.0-gamma)/rho );
        D2flux.set( i, 2, 3, 1, gm1*u2/rho );    D2flux.set( i, 2, 1, 3, gm1*u2/rho );
        D2flux.set( i, 2, 3, 3, -gm1/rho );
        D2flux.set( i, 2, 4, 1, gm1*u3/rho ); D2flux.set( i, 2, 1, 4, gm1*u3/rho );
        D2flux.set( i, 2, 4, 4, -gm1 / rho ); 

        // Derivatives with respect to f3:
        D2flux.set( i, 3, 1, 1, 2.0*u1*u2/rho );
        D2flux.set( i, 3, 2, 1, -u2/rho       );  D2flux.set( i, 3, 1, 2, -u2/rho     );  
        D2flux.set( i, 3, 3, 1, -u1/rho       );  D2flux.set( i, 3, 1, 3, -u1/rho     );
        D2flux.set( i, 3, 3, 2, 1.0/rho       );  D2flux.set( i, 3, 2, 3, 1.0/rho     );  

        
        // Derivatives with respect to f4:
        D2flux.set( i, 4, 1, 1, 2.0*u1*u3/rho );
        D2flux.set( i, 4, 2, 1, -u3/rho       ); D2flux.set( i, 4, 1, 2, -u3/rho       );
        D2flux.set( i, 4, 4, 1, -u1 / rho     ); D2flux.set( i, 4, 1, 4, -u1 / rho     );
        D2flux.set( i, 4, 4, 2, 1.0/rho       ); D2flux.set( i, 4, 2, 4, 1.0/rho       );
         


        // Derivatives with respect to f5:
        D2flux.set(i,5,1,1, u1*( 2.0*gamma*energy - 3.0*rho*gm1*umag2 ) / (rho*rho) );

        // TODO - this one is NOT correct:
        D2flux.set(i,5,2,1, -( gamma*energy/rho - (gm1)*(umag2) - 2.0*gm1*u1*u1 )/rho );
        D2flux.set(i,5,1,2, D2flux.get(i,5,2,1) );

        D2flux.set(i,5,2,2, -3.0*gm1*u1/rho );

        D2flux.set(i,5,3,1, 2.0*u1*gm1*u2/rho );
        D2flux.set(i,5,1,3, D2flux.get(i,5,3,1) );

        D2flux.set(i,5,3,2, -gm1*u2/rho );
        D2flux.set(i,5,2,3, D2flux.get(i,5,3,2) );

        D2flux.set(i,5,3,3, -gm1*u1/rho );

        D2flux.set(i,5,4,1, 2.0*u1*gm1*u3/rho );
        D2flux.set(i,5,1,4, D2flux.get(i,5,4,1) );

        D2flux.set(i,5,4,2, -u3*gm1/rho       );
        D2flux.set(i,5,2,4, D2flux.get(i,5,4,2) );

        D2flux.set(i,5,4,4, -gm1*u1/rho       );

        D2flux.set(i,5,5,1, -u1*gamma/rho     );
        D2flux.set(i,5,1,5, D2flux.get(i,5,5,1) );

        D2flux.set(i,5,5,2, gamma/rho         );
        D2flux.set(i,5,2,5, D2flux.get(i,5,5,2) );

    }

}
