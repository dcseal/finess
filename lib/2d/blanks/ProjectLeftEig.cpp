#include "tensors.h"

// * TEMPLATE *
//
// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
// For the 2D balance law: 
//
//            q_t + f_x + g_y = psi
//
// There are two potential matrices to consider: f'(q) and g'(q).
//
//   ixy == 1 refers to the eigenstructure of f'(q), and
//   ixy == 2 refers to the eigenstructure of g'(q).
//
// See also: ProjectRightEig.
void ProjectLeftEig(int ixy, 
    const dTensor1& Aux_ave, const dTensor1& Q_ave, const dTensor2& Qvals,
		  dTensor2& Wvals)
{

  const int meqn = Qvals.getsize(1);
  const int mpts = Qvals.getsize(2);

}
