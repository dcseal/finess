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
 *     a.) Call FinSolve[TS-method], where TS-method is a valid time-stepping
 *     option.  (e.g. FinSolveRK, FinSolveSDC, FinSolveUser).
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
#include "StateVars.h"
#include "RunFinpack.h"           // Function declarations

int RunFinpack( )
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
    const int&     meqn     = global_ini_params.get_meqn();
    const int&     maux     = global_ini_params.get_maux();
    const int&     mx       = global_ini_params.get_mx();
    const int&     mbc      = global_ini_params.get_mbc();
    const double&  xlow     = global_ini_params.get_xlow();
    const double&  xhigh    = global_ini_params.get_xhigh();
    const double&  dx       = global_ini_params.get_dx();

    // Dimension arrays
    StateVars Qstate(0., mx, meqn, maux, mbc );

    // Set any auxiliary variables on computational grid
    // Set values and apply L2-projection
    if(maux >0)
    {  SampleFunction(1-mbc, mx+mbc, Qstate.ref_q(), Qstate.ref_aux(), Qstate.ref_aux(), &AuxFunc);  }

    // Set initial data on computational grid
    // Set values and apply L2-projection
    // SampleFunction( 1-mbc, mx+mbc, qnew, aux, qnew, &QinitFunc);
    SampleFunction( 1-mbc, mx+mbc, Qstate.ref_q(), Qstate.ref_aux(), Qstate.ref_q(), &QinitFunc);

    // Run AfterQinit to set any necessary variables
    AfterQinit( Qstate );

    // Output initial data to file
    // For each element, we output ``method[1]'' number of values
    Output( Qstate, 0 );

    // Compute conservation and print to file
    ConSoln( Qstate );

    // Main loop for time stepping
    double tstart = 0.0;
    double tend   = 0.0;
    double dtout  = tfinal/double(nout);    
    for (int n=1; n<=nout; n++)
    {        
        tstart = tend;      assert_lt( fabs(Qstate.get_t()-tend), 1e-13 );
        tend   = tstart + dtout;

        // Solve hyperbolic system from tstart to tend
        if (time_stepping_method == IniParams::TimeSteppingMethod::RK)
        {  
            // Runge-Kutta time-stepping scheme
            FinSolveRK( Qstate, tend, dtv );
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::LxW)
        {
            // Lax-Wendroff time stepping
            FinSolveLxW( Qstate, tend, dtv );
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::MD)
        { 
            // Multiderivative time-stepping
            FinSolveMD( Qstate, tend, dtv ); 
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::SDC)
        { 
            // Spectral deferred correction time-stepping
            FinSolveSDC( Qstate, tend, dtv ); 
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::USER_DEFINED)
        {
            // User-defined time-stepping scheme
            FinSolveUser( Qstate, tend, dtv );
        }
        else
        {
            printf("Error: RunFinpack.  Time stepping method not implemented\n" );
            exit(1);
        }

        // Output data to file
        Output( Qstate, n );

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
    QinitFunc(xpts, qvals);
}

void AuxFunc(const dTensor1& xpts, const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, dTensor2& auxvals)
{
    void AuxFunc(const dTensor1& xpts, dTensor2& auxvals);
    AuxFunc(xpts, auxvals);
}
