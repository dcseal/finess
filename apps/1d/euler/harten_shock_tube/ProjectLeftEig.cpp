#include <cmath>
#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
// Input:
//
// Aux_ave( 1:maux )
//   Q_ave( 1:meqn )
//   Qvals( 1:meqn, 1:numpts )
//
// Output:
//
//   Wvals( 1:meqn, 1:numpts )
//

///@brief User-supplied that transforms physical quantities into characteristic variables
///
///@param Aux_ave  One-dimensional array.  Index range: {1, ..., maux}.
///                Stores @f$ \mathrm{aux}_{i-1/2} @f$, the average of aux[i, :] and aux[i-1, :]
///@param Q_ave    One-dimensional array.  Index range: {1, ..., meqn}.
///                Stores @f$ q_{i-1/2} @f$, the average of q[1, :] and q[i-1, :]
///@param Qvals    Two-dimensional array.  Index range: {1, ..., meqn}*{1, ..., 2r}
///                Physical quantities (q, or f) at @f$ x_{i-r}, ..., x_{i+(r-1)}@f$.
///@param Wvals    Two-dimensional array to store the characteristic variables.
///                Index range: {1, ..., meqn}*{1, ..., 2r}
///                Will hold w-values (if Qvals holds q-values),
///                                 or g-values (ifQvals hold f-values).
///
///Conceptually, does the following:
///-# Compute @f$ f'(q_{i-1/2}) @f$ based on <tt>Aux_ave</tt> and <tt>Q_ave</tt>.
///-# Decompose @f$ f'(q_{i-1/2}) = R\Lambda R^{-1}@f$, where @f$ \Lambda @f$ is diagonal.
///-# @f$ \textrm{Wvals} \leftarrow R^{-1} \cdot \textrm{Qvals} @f$
///
///Implementation is more compact than the conceptual description.
void ProjectLeftEig(
    const dTensor1& Aux_ave, 
    const dTensor1& Q_ave,
    const dTensor2& Qvals,
    dTensor2& Wvals)
{    

    const int meqn = Qvals.getsize(1);
    const int npts = Qvals.getsize(2);

    // Average states
    const double gamma  = eulerParams.gamma;
    const double rho    = Q_ave.get(1);
    const double u1     = Q_ave.get(2)/rho;
    const double u2     = Q_ave.get(3)/rho;
    const double u3     = Q_ave.get(4)/rho;
    const double energy = Q_ave.get(5);
    const double umag2  = (u1*u1 + u2*u2 + u3*u3);
    const double press  = (gamma-1.0e0)*(energy-0.5e0*rho*umag2);
    const double c      = sqrt(fabs(gamma*press/rho));
    const double H      = (energy+press)/rho; 

    // Project onto left eigenvectors
    for (int k=1; k <= npts; k++)
    {
        Wvals.set(1,k, ((umag2/2.0-H-u1*c)*Qvals.get(2,k) 
                    + (umag2/2.0*(c-u1)+H*u1)*Qvals.get(1,k) 
                    + c*(Qvals.get(5,k)-u2*Qvals.get(3,k)
                        - u3*Qvals.get(4,k)))/(c*(2.0*H-umag2)) );

        Wvals.set(2,k, 2.0*((H-umag2)*Qvals.get(1,k) + u1*Qvals.get(2,k) 
                    + u2*Qvals.get(3,k) + u3*Qvals.get(4,k) 
                    - Qvals.get(5,k))/(2.0*H-umag2) );

        Wvals.set(3,k, Qvals.get(3,k)-u2*Qvals.get(1,k) );

        Wvals.set(4,k, Qvals.get(4,k)-u3*Qvals.get(1,k) );

        Wvals.set(5,k, ((H-umag2/2.0-u1*c)*Qvals.get(2,k) 
                    + (umag2/2.0*(c+u1)-H*u1)*Qvals.get(1,k)
                    + c*(Qvals.get(5,k)-u2*Qvals.get(3,k) 
                        - u3*Qvals.get(4,k)))/(c*(2.0*H-umag2)) );
    }
}
