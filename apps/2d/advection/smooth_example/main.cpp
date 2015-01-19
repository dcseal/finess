// =========================================================================
//
//  --------------------------------------------
//  FINESS: See LICENSE.txt for licensing details
//  --------------------------------------------
//
//    Lead Developer:  
//             David C. Seal
//             Michigan University
//             Department of Mathematics
//             619 Red Cedar Road
//             East Lansing, MI 48823
//             seal@math.msu.edu
//
// =========================================================================

int main(int argc, char* argv[])
{

    // NOTE: You should not have to modify this part of the code.
    //
    //       To change parameters, modify the following files:
    //            1. parameters.ini -- basic data file, can modify
    //                  # of grid points, time step, order of
    //                  accuracy in both space and time, etc...
    //            2. QinitFunc.cpp -- initial condition file
    //            3. AuxFunc.cpp -- auxiliary variable file
    //            4. SourceTermFunc.pp -- source term file
    //            5. FluxFunc.cpp -- flux function file
    //            6. SetWaveSpd.cpp -- eigenvalues of flux Jacobian file
    //            7. ProjectLeftEig.cpp -- left eigenvectors of flux Jacobian file
    //            8. ProjectLeftEig.cpp -- right eigenvectors of flux Jacobian file
    //            9. SetBndValues.cpp -- boundary conditions files
    //

    int main_global(int argc, char* argv[]);
    return main_global(argc,argv);

}
