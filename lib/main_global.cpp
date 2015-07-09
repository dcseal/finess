#include <string>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "dogdefs.h"
#include "dimdefs.h"
#include "IniParams.h"


#include <fenv.h>

/*
 * Common main function that's called by every application.
 *
 * Each application has its own main file, and Makefile that builds that local
 * main.  When the code is built in any application directory, every
 * application links to this common main function, that executes the code.
 *
 * The purpose of placing this extra layer between appname/main.cpp and
 * RunFinpack is to make the main function in each application as short as
 * possible.
 *
 */
int main_global(int argc, char* argv[])
{
//   feenableexcept(FE_INVALID  |
//                FE_DIVBYZERO |
//                FE_OVERFLOW  );
//     feenableexcept(FE_INVALID  |
//                  FE_DIVBYZERO |
//                  FE_OVERFLOW  |
//                  FE_UNDERFLOW);
    
    using std::string;
    using std::cout;
    using std::setprecision;
    using std::setw;
    using std::scientific;
    using std::endl;

    string parameters_ini_filename;

    parameters_ini_filename = 
        argc == 1 ? "parameters.ini" : argv[1];
       
    cout << "Running with configuration file: "
         << parameters_ini_filename
         << endl;

    global_ini_params.init(parameters_ini_filename);

    // Get current time
    double time1 = time(NULL);

    // Call startscript (Default: scripts/startscript, otherwise it checks for
    // a local file called 'startscript' from the application's directory)
    void RunStartScript(string parameters_ini_filename, string finess_exe_path);
    RunStartScript(parameters_ini_filename, argv[0]);

    // Call the ``RunFinpack'' routine, which executes the code
    // Each dimension has its own version of this routine.
    int RunFinpack(string parameters_ini_filename );
    int m = RunFinpack(parameters_ini_filename);

    // Get current time
    double time2 = time(NULL);

    // Output elapsed running time
    cout << setprecision(5);
    cout << " Total elapsed time in seconds = " << setw(11) 
        << scientific << time2-time1 << endl << endl;

    return m;

}
