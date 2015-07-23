
// This is a user-supplied routine that projects
// W onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Q
//
class dTensor1;
class dTensor2;
void ProjectRightEig( const dTensor1& Aux_ave,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{    

    // We do not use this in the implicit solver
    // retained to ensure compatibility with general FINESS structure

}
