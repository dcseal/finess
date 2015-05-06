#include "Components.h"
// This is a user-supplied routine that projects
// Q onto the left eigenvectors of the flux 
// Jacobian; the result is stored in W
//
class dTensor1;
class dTensor2;
void ProjectLeftEig( int ixy, const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{    
    void ProjectLeftEig_FiveMoment( int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
    ProjectLeftEig_FiveMoment(ixy, 0, Q_ave, Q, W);
    ProjectLeftEig_FiveMoment(ixy, 5, Q_ave, Q, W);

    void ProjectLeftEig_Maxwell(int ixy, int n_offset,
        const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
    ProjectLeftEig_Maxwell(ixy, 10, Q_ave, Q, W);

//  TODO - I'm not sure what this part is for. -DS
//  void ProjectLeftConvectedScalar(int idx,
//      const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W);
//  ProjectLeftConvectedScalar(_entropy_i, Q_ave, Q, W);
//  ProjectLeftConvectedScalar(_entropy_e, Q_ave, Q, W);

}
