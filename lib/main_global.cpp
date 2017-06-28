#include <string>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "dogdefs.h"
#include "dimdefs.h"
#include "IniParams.h"

/*
 * Common main function that's called by every application.
 *
 * Each application has its own main and Makefile that builds a local
 * main function.  When the code is built in any application directory, every
 * application links to this common main function, which serves as the entrance
 * to being able to execute the FINESS software.
 *
 * The purpose of placing this extra layer between $(FINESS)/[appname]/main.cpp and
 * RunFinpack is to make the main function in each application as short as
 * possible.
 *
 * Note that the default parameters file is parameters.ini but additional
 * options can be given on the command line via:
 *
 *      ./finess [parameters_file.ini]
 *
 * This is a helpful way to, for example, perform a batch script parameter
 * study. (Just make sure to change the output_dir parameter in the
 * parameters.ini file that you link to.)
 *
 */
int main_global(int argc, char* argv[])
{

    // IO routines to print helpful messages to the user
    using std::string;
    using std::cout;
    using std::setprecision;
    using std::setw;
    using std::scientific;
    using std::endl;

    string parameters_ini_filename;

    // Choose what parameters file to run the code with.  The default setting
    // is to use 'parameters.ini.'  If another argument is given in the command
    // line, then that parameters file is run instead.
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
    void RunStartScript(string parameters_ini_filename);
    RunStartScript(parameters_ini_filename);

    // Call the ``RunFinpack'' routine, which executes the code
    // Each dimension has its own version of this routine.
    int RunFinpack( );
    int m = RunFinpack( );

    // Get current time
    double time2 = time(NULL);

    // Output elapsed running time
    cout << setprecision(5);
    cout << " Total elapsed time in seconds = " << setw(11) 
        << scientific << time2-time1 << endl << endl;

    return m;

}
