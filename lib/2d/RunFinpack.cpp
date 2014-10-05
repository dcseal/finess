#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "dogdefs.h"
#include "IniParams.h"            // accessors for the parameters.ini file
#include "RunFinpack.h"           // Function declarations
#include "StateVars.h"


/*
 * Top level function to RunFinpack.  Briefly, this function calls the
 * following functions in the following order:
 *
 * TODO
 *
 */
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
    const int&     my       = global_ini_params.get_my();
    const int&     mbc      = global_ini_params.get_mbc();

    // Dimension arrays
    StateVars Qnew(0., mx, my, meqn, maux, mbc );
    dTensorBC3& qnew = Qnew.ref_q();
    dTensorBC3& aux  = Qnew.ref_aux();

    // Set any auxiliary variables on computational grid
    // Set values and apply L2-projection
    if( maux > 0 )
    {  SampleFunction(1-mbc, mx+mbc, 1-mbc, my+mbc, qnew, aux, aux, &AuxFunc);  }

    // Set initial data on computational grid
    // Set values and apply L2-projection
    SampleFunction( 1-mbc, mx+mbc, 1-mbc, my+mbc, qnew, aux, qnew, &QinitFunc);

    // Run AfterQinit to set any necessary variables
    AfterQinit( Qnew );

    // Output initial data to file
    // For each element, we output ``method[1]'' number of values
    Output( Qnew, 0 );

    // Compute conservation and print to file
    ConSoln( Qnew );

    // Main loop for time stepping
    double tstart = 0.0;
    double tend   = 0.0;
    double dtout  = tfinal/double(nout);    
    for (int n=1; n<=nout; n++)
    {        
        tstart = tend;      
        tend   = tstart + dtout;

        // Solve hyperbolic system from tstart to tend
        if (time_stepping_method == IniParams::TimeSteppingMethod::RK)
        {  
            // Runge-Kutta time-stepping scheme
            FinSolveRK( Qnew, tend, dtv );
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::LxW)
        {
            // Lax-Wendroff time stepping
            FinSolveLxW(Qnew, tend, dtv );
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::MD)
        {
            // Multiderivative time stepping
            FinSolveMD(Qnew, tend, dtv );
        }
        else if (time_stepping_method == IniParams::TimeSteppingMethod::USER_DEFINED)
        {
            // User-defined time-stepping scheme
            FinSolveUser(Qnew, tend, dtv );
        }
        else
        {
            printf("Time stepping method not implemented\n");
            exit(1);
        }

        // Output data to file
        Output( Qnew, n );

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
