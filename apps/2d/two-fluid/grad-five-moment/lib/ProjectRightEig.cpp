#include "Components.h"
// This is a user-supplied routine that projects
// W onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Q
//
class dTensor1;
class dTensor2;
void ProjectRightEig( int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{    
    void ProjectRightEig_FiveMoment(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
    ProjectRightEig_FiveMoment(ixy, 0, Q_ave, W, Q);
    ProjectRightEig_FiveMoment(ixy, 5, Q_ave, W, Q);

    void ProjectRightEig_Maxwell(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
    ProjectRightEig_Maxwell(ixy, 10, Q_ave, W, Q);

//  TODO - I'm not sure what this part is for.  -DS
//  void ProjectRightConvectedScalar(int idx,
//      const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q);
//  ProjectRightConvectedScalar(_entropy_i, Q_ave, W, Q);
//  ProjectRightConvectedScalar(_entropy_e, Q_ave, W, Q);

}
