#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "dogdefs.h"
#include "dimdefs.h"


#include "IniParams.h"
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
    using std::string;
    using std::cout;
    using std::setprecision;
    using std::setw;
    using std::scientific;
    using std::endl;

    global_ini_params.init("parameters.ini");


    // Get current time
    double time1 = time(NULL);


    // Call startscript (Default: scripts/startscript, otherwise it checks for
    // a local file called 'startscript' from the application's directory)
    void RunStartScript(int ndims);
    RunStartScript(NDIMS);

    // Call the ``RunFinpack'' routine, which executes the code
    // Each dimension has its own version of this routine.
    int RunFinpack( string outputdir );
    int m = RunFinpack( global_ini_params.get_output_dir() );

    // Get current time
    double time2 = time(NULL);

    // Output elapsed running time
    cout << setprecision(5);
    cout << " Total elapsed time in seconds = " << setw(11) 
        << scientific << time2-time1 << endl << endl;

    return m;

}
