#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "dogdefs.h"
#include "IniDocument.h"
#include "DogSolver.h"

/*
 * Common main function that's called by every (1D) application.
 *
 * Each application has its own main file, and Makefile that builds that local
 * main.  When the code is built in any application directory, every (1D) 
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

    DogSolver::parse_arguments(argc,argv);

    void RunStartScript(int ndims);
    RunStartScript(1);

    // Call the ``RunFinpack'' routine, which executes the code
    int RunFinpack(string outputdir);
    int m = RunFinpack(get_outputdir());

    // Get current time
    double time2 = time(NULL);

    // Output elapsed running time
    cout << setprecision(5);
    cout << " Total elapsed time in seconds = " << setw(11) 
        << scientific << time2-time1 << endl << endl;

    return m;
}
