/*
 * Top level function to RunFinpack.  Briefly, this function calls the
 * following functions in the following order:
 *
 * 1.)  Iniitalize global structs dogParams and dogParamsCart1
 *
 * 2.) Call InitApp.  (additional application specific parameters)
 *
 * 3.) Write qhelp.dat to the output directory.  This is a total of two
 * functions calls: one on global_ini_params.write_qhelp and one on
 * global_ini_params.write_qhelp.
 *
 * 4.) Call L2Project.  This Projects initial conditions onto basis functions.
 *
 * 5.) Call AfterQinit.  This is called once per simulation and can be used to
 * set up extra variables.
 *
 * 6.) Call Output - output initial conditions to the output directory.
 *
 * 7.) Call ConSoln - this is a call for saving 'conserved' quantities.  This
 * function is called once per time step.
 *
 * 8.) Run the main time stepping loop.  This consists of calling the
 * following two functions, once for each frame the user requested:
 *
 *     a.) Call DogSolve[TS-method], where TS-method is a valid time-stepping
 *     option.  (e.g. DogSolveRK, DogSolveSDC, DogSolveUser).
 *
 *     b.) Call Output to print data to file
 *
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "dogdefs.h"
#include "IniParams.h"            // accessors for the parameters.ini file
#include "RunFinpack.h"           // Function declarations

int RunFinpack(std::string outputdir)
{

    using std::cout;
    using std::endl;
    using std::string;
    using std::scientific;
    using std::setw;
    using std::setprecision;

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

    // Get parameters and print to screen
    cout << global_ini_params.ini_doc_as_string() << endl;
    const IniParams::TimeSteppingMethod::enum_type time_stepping_method = 
	  global_ini_params.get_time_stepping_method();

    const int&     nout     = global_ini_params.get_nout();
    const double&  tfinal   = global_ini_params.get_tfinal();
    double dtv[2+1];
    dtv[1] = global_ini_params.get_initial_dt();
    dtv[2] = global_ini_params.get_max_dt();
    const double nonsense_double = 0;
    const double   cflv[]   = {nonsense_double, global_ini_params.get_max_cfl(), global_ini_params.get_desired_cfl()};
    const int      nv       = global_ini_params.get_nv();
    const int&     meqn     = global_ini_params.get_meqn();
    const int&     maux     = global_ini_params.get_maux();
    const int&     mdim     = global_ini_params.get_ndims();     assert_eq( mdim, 1 );
    const int&     mx       = global_ini_params.get_mx();
    const int&     mbc      = global_ini_params.get_mbc();
    const double&  xlow     = global_ini_params.get_xlow();
    const double&  xhigh    = global_ini_params.get_xhigh();
    const double&  dx       = global_ini_params.get_dx();
//  const int&     mrestart = global_ini_params.get_mrestart();

    // Output helpful stuff to qhelp.dat for plotting purposes
    // This is deprecated (replaced with global IniParams)
//  string qhelp;
//  qhelp = outputdir+"/qhelp.dat";
//  global_ini_params.write_qhelp(qhelp.c_str());
//  global_ini_params.append_qhelp(qhelp.c_str());

    // Dimension arrays
    dTensorBC2    qnew(mx, meqn, mbc);
    dTensorBC2    qold(mx, meqn, mbc);
    dTensorBC1    smax(mx, mbc);
    dTensorBC2    aux (mx, iMax(maux, 1), mbc);

    // Set any auxiliary variables on computational grid
    // Set values and apply L2-projection
    if(maux >0)
    {  SampleFunction(1-mbc, mx+mbc, qnew, aux, aux, &AuxFunc);  }

    // Set initial data on computational grid
    // Set values and apply L2-projection
    SampleFunction( 1-mbc, mx+mbc, qnew, aux, qnew, &QinitFunc);

    // Run AfterQinit to set any necessary variables
    AfterQinit( aux, qnew);

    // Output initial data to file
    // For each element, we output ``method[1]'' number of values
    Output( aux, qnew, 0.0, 0, outputdir);

    // Compute conservation and print to file
    ConSoln( aux, qnew, 0.0, outputdir);

    // Main loop for time stepping
    double tstart = 0.0;
    double tend   = 0.0;
    double dtout = tfinal/double(nout);    
    for (int n=1; n<=nout; n++)
    {        
        tstart = tend;      
        tend = tstart + dtout;

        // Solve hyperbolic system from tstart to tend
        if (time_stepping_method == IniParams::TimeSteppingMethod::RK)
        {  
            // Runge-Kutta time-stepping scheme
            FinSolveRK( aux, qold, qnew, smax, tstart, tend, 
                    nv, dtv, cflv, outputdir);
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::LxW)
        {
            FinSolveLxW( 
                aux, qold, qnew, smax, tstart, tend, 
                nv, dtv, cflv, outputdir);
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::MD)
        {
            FinSolveMD( 
                aux, qold, qnew, smax, tstart, tend, 
                nv, dtv, cflv, outputdir);
        }
//      else if( time_stepping_method == "SDC")
//      {
//          FinSolveSDC( 
//              aux, qold, qnew, smax, tstart, tend, 
//              nv, dtv, cflv, outputdir);
//      }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::USER_DEFINED)
        {
            // User-defined time-stepping scheme
            DogSolveUser( aux, qold, qnew, smax, tstart, tend, 
                    nv, dtv, cflv, outputdir);
        }
        else
        {
            printf("Error: RunFinpack.  Time stepping method not implemented\n" );
            exit(1);
        }

        // Output data to file
        Output( aux, qnew, tend, n, outputdir);

        // Done with solution from tstart to tend
        cout << setprecision(5);
        cout << "FINESS: Frame " << setw(3) << n;
        cout << ": plot files done at time t =";
        cout << setw(12) << scientific << tend << endl;
        cout << endl;
    }

    return 0;
}

// Wrapper functions to make the calls to Qinit and AuxFunc make sense when
// passed into SampleFunction
void QinitFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& qvals)
{
    void QinitFunc(const dTensor1& xpts, dTensor2& qvals);
    QinitFunc(xpts,qvals);
}

void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals)
{
    void AuxFunc(const dTensor1& xpts, dTensor2& auxvals);
    AuxFunc(xpts,auxvals);
}
