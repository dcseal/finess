/*
 * Top level function to RunFinpack.  Briefly, this function calls the
 * following functions in the following order:
 *
 * TODO
 *
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"            // accessors for the parameters.ini file
#include "DogParamsCart2.h"       // accessors for the parameters.ini file
#include "IniDocument.h"
#include "RunFinpack.h"           // Function declarations

int RunFinpack(string outputdir)
{

    // Output title information
    cout << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << "   | FINESS:  The FINite difference ESSentially   |   " << endl;
    cout << "   |          non-oscillatory software Package    |   " << endl;
    cout << "   | Developed by the research group of           |   " << endl;
    cout << "   |            David C. Seal                     |   " << endl;
    cout << "   |            Department of Mathematics         |   " << endl;
    cout << "   |            Michigan State University         |   " << endl;
    cout << "   ------------------------------------------------   " << endl;
    cout << endl;

    // Get parameters
    dogParams.init();
    dogParamsCart2.init(ini_doc);
    cout << endl;

    // Get addtional parameters
//  InitApp( ini_doc );
//  cout << endl;

    const string time_stepping_method = dogParams.get_time_stepping_method();
    const int&     nout     = dogParams.get_nout();
    const double&  tfinal   = dogParams.get_tfinal();
    double dtv[2+1];
    dtv[1] = dogParams.get_initial_dt();
    dtv[2] = dogParams.get_max_dt();
    const double*  cflv     = dogParams.get_cflv();
    const int      nv       = dogParams.get_nv();
    const int&     meqn     = dogParams.get_meqn();
    const int&     maux     = dogParams.get_maux();
    const int&     mdim     = dogParams.get_ndims();     assert_eq( mdim, 2 );
    const int&     mx       = dogParamsCart2.get_mx();
    const int&     my       = dogParamsCart2.get_my();
    const int&     mbc      = dogParamsCart2.get_mbc();
    const int&     mrestart = dogParams.get_mrestart();
    const int        mnodes = mx + 1;

    // Output helpful stuff to qhelp.dat for plotting purposes
    string qhelp;
    qhelp=outputdir+"/qhelp.dat";
    dogParams.write_qhelp(qhelp.c_str());
    dogParamsCart2.append_qhelp(qhelp.c_str());


    // Dimension arrays
    dTensor3      node(mx, my, mdim);
    dTensor2      prim_vol(mx, my);
    dTensorBC3    qnew(mx, my, meqn, mbc);
    dTensorBC3    qold(mx, my, meqn, mbc);
    dTensorBC2    smax(mx, my, mbc);
    dTensorBC3    aux (mx, my, iMax(maux, 1), mbc);

    // Construct 1D grid (trivial for uniform case)
    GridSetup( node, prim_vol);

/*
    // Set any auxiliary variables on computational grid
    // Set values and apply L2-projection
    if(maux >0)
    {  SampleFunction(1-mbc, mx+mbc, node, qnew, aux, aux, &AuxFunc);  }

    // Set initial data on computational grid
    // Set values and apply L2-projection
    SampleFunction( 1-mbc, mx+mbc, node, qnew, aux, qnew, &QinitFunc);

    // Run AfterQinit to set any necessary variables
    AfterQinit(node, aux, qnew);

    // Output initial data to file
    // For each element, we output ``method[1]'' number of values
    Output(node, aux, qnew, 0.0, 0, outputdir);

    // Compute conservation and print to file
    ConSoln( node, aux, qnew, 0.0, outputdir);

    // Main loop for time stepping
    double tstart = 0.0;
    double tend   = 0.0;
    double dtout = tfinal/double(nout);    
    for (int n=1; n<=nout; n++)
    {        
        tstart = tend;      
        tend = tstart + dtout;

        // Solve hyperbolic system from tstart to tend
        if (time_stepping_method == "Runge-Kutta")
        {  
            // Runge-Kutta time-stepping scheme
            FinSolveRK(node, prim_vol, aux, qold, qnew, smax, tstart, tend, 
                    nv, dtv, cflv, outputdir);
        }
        else if (time_stepping_method == "User-Defined")
        {
            // User-defined time-stepping scheme
            DogSolveUser(node,  prim_vol, aux, qold, qnew, smax, tstart, tend, 
                    nv, dtv, cflv, outputdir);
        }

        // Output data to file
        Output(node, aux, qnew, tend, n, outputdir);

        // Done with solution from tstart to tend
        cout << setprecision(5);
        cout << "FINESS: Frame " << setw(3) << n;
        cout << ": plot files done at time t =";
        cout << setw(12) << scientific << tend << endl;
        cout << endl;
    }
*/
    return 0;
}

// Wrapper functions to make the calls to Qinit and AuxFunc make sense when
// passed into SampleFunction
void QinitFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals)
{
    void QinitFunc(const dTensor2& xpts, dTensor2& qvals);
    QinitFunc(xpts,qvals);
}

void AuxFunc(const dTensor2& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals)
{
    void AuxFunc(const dTensor2& xpts, dTensor2& auxvals);
    AuxFunc(xpts, auxvals);
}
