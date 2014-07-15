#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "dogdefs.h"
#include "IniDocument.h"
#include "OutputDir.h"
#include "dimdefs.h"

/*
 * Common main function that's called by every (3D) application.
 *
 * Each application has its own main file, and Makefile that builds that local
 * main.  When the code is built in any application directory, every (3D) 
 * application links to this common main function, that executes the code.
 *
 * The purpose of placing this extra layer between appname/main.cpp and
 * RunFinpack is to make the main function in each application as short as
 * possible.
 *
 */
int main_global(int argc, char* argv[])
{
    // Open parameters.ini file, (and read [dogParams]? -DS)
    ini_doc.initFromFile("parameters.ini");
    IniDocument::Section& ini_sec = ini_doc["dogParams"];

    // Get current time
    double time1 = time(NULL);

    // Parse the command line arguments.  (e.g. -o sets a different output
    // directory, and -d sets a different debug level).
    //
    // TODO - there is no need for this to be a class.  A single global
    // variable containing the output directory would suffice. (-DS)
    OutputDir::parse_arguments(argc, argv);

    // Call startscript (Default: scripts/startscript, otherwise it checks for
    // a local file called 'startscript' from the application's directory)
    void RunStartScript(int ndims);
    RunStartScript(NDIMS);

    // Call the ``RunFinpack'' routine, which executes the code
    int RunFinpack( string outputdir );
    int m = RunFinpack( get_outputdir() );

    // Get current time
    double time2 = time(NULL);

    // Output elapsed running time
    cout << setprecision(5);
    cout << " Total elapsed time in seconds = " << setw(11) 
        << scientific << time2-time1 << endl << endl;

    return m;
}
