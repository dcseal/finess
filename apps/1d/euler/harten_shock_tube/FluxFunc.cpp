#include "dogdefs.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Euler equations for gas dynamics
//
///@brief Sets flux function for Euler equations.
///
///- Euler equations for gas dynamics:
/**
@f[
\left(\begin{array}{c}
\rho\\
\rho u_{1}\\
\rho u_{2}\\
\rho u_{3}\\
E
\end{array}\right)_{t}+\left(\begin{array}{c}
\rho u_{1}\\
\rho u_{1}+p\\
\rho u_{1}u_{2}\\
\rho u_{1}u_{3}\\
(E+p)u_{1}
\end{array}\right)_{x}=0@f]
*/
///   where @f$ p = (\gamma-1) (E - \frac{1}{2} \rho (u_1^2 + u_2^2 + u_3^2))  @f$
///
///- See #SampleFunction(...)
///  for possible calls through function pointer to the current function.
///
///@param xpts Supposed to be the x-coordinates of mesh points.  Not used in this problem.
///@param Q Two-dimensional array with dimension 2 of size 5.  Represents the state vector of solution,
/** i.e. @f$ \left(\begin{array}{c}
\rho\\
\rho u_{1}\\
\rho u_{2}\\
\rho u_{3}\\
E
\end{array}\right)@f$  at mesh points.*/
/// For example, Q[i, 1] is @f$\rho@f$ value at i-th mesh point, etc.
///@param Aux Two-dimensional array containing "auxiliary" values.  Not used in this problem.
///@param flux Sink for the "state" vector for flux function, i.e. where to put values of 
/** 
@f$
\left(\begin{array}{c}
\rho u_{1}\\
\rho u_{1}+p\\
\rho u_{1}u_{2}\\
\rho u_{1}u_{3}\\
(E+p)u_{1}
\end{array}\right)
@f$  (which are evaluated at mesh points)
*/
///
///
///@todo Clarify use of <tt>xpts</tt>, and <tt>Aux</tt>. --- Why are they necessary?
///@todo Clarify the form of the equations 
///      --- this is different from, e.g. LeVeque Finite Volume Methods Book, page 293, (14.8)
void FluxFunc(const dTensor1& xpts, 
	      const dTensor2& Q,
	      const dTensor2& Aux, 
	      dTensor2& flux)
{
  const int numpts=xpts.getsize();

  // Gas constant
  const double gamma = eulerParams.gamma;

  for (int i=1; i<=numpts; i++)
    {
      double x = xpts.get(i);      
    
      // Variables
      const double rho    = Q.get(i,1);
      const double u1     = Q.get(i,2)/rho;
      const double u2     = Q.get(i,3)/rho;
      const double u3     = Q.get(i,4)/rho;
      const double energy = Q.get(i,5);
      const double press  = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));
    
      // Flux function
      flux.set(i,1, rho*u1 );
      flux.set(i,2, rho*u1*u1 + press );
      flux.set(i,3, rho*u1*u2 );
      flux.set(i,4, rho*u1*u3 );
      flux.set(i,5, u1*(energy+press) ); 
    }
}
