#include "tensors.h"

// *REQUIRED*
//
// This is a user-supplied routine sets the maximum wave speeds for use in a
// Riemann solver.
//
// For the 3D balance law: 
//
//            q_t + f_x + g_y + h_z = psi, or rather, q_t + div( R ) = 0.
//
// Given a normal vector n, the matrix n \cdot R is diagonalizable into
//
//          n \cdot R = R \Lambda R^{-1}.
//
// This routine computes s1 = min( \Lambda ) and s2 = max( \Lambda ).
//
// Input:
//
//    nvec( 1:ndim )  - Normal vector
//   xedge( 1:ndim )  - Location of the Riemann problem
//      Ql( 1:meqn )  - Left state
//      Qr( 1:meqn )  - Right state
//    Auxl( 1:maux )  - auxiliary function evaluated at the left state
//    Auxr( 1:maux )  - auxiliary function evaluated at the right state
//
// Output:
//
//   s1 = Min( \Lambda )
//   s2 = Max( \Lambda )
//
// See also: ProjectLeftEig, ProjectRightEig. 
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        double& s1,double& s2)
{ }
