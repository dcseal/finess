#include "tensors.h"

// * TEMPLATE *
//
// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals.
//
// For the 3D balance law: 
//
//            q_t + f_x + g_y + h_z = psi
//
// There are three potential matrices to consider: f'(q), g'(q) and h'(q)
//
//   ixy == 1 refers to the eigenstructure of f'(q), and
//   ixy == 2 refers to the eigenstructure of g'(q).
//   ixy == 3 refers to the eigenstructure of h'(q).
//
// See also: ProjectLeftEig.
void ProjectRightEig(int ixy, 
    const dTensor1& Aux_ave, const dTensor1& Q_ave, const dTensor2& Wvals,
    dTensor2& Qvals)
{ 

  const int meqn = Qvals.getsize(1);
  const int mpts = Qvals.getsize(2);

}
