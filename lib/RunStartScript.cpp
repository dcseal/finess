#include <string>
#include <stdlib.h> // for system()
#include "assert.h"
#include "debug.h"

#include "IniParams.h"

// Function used to call startscript from $(FINESS)/scripts;
//
// That shell script creates the output directory and copies a few files in
// the output folder
void RunStartScript(int ndims )
{
    char command_str[1024];
//    const char* get_outputdir();
    std::string output_dir = global_ini_params.get_output_dir();

    // run startscript
    // to create output directory if it does not exist
    // and copy data files to output directory
    int numchars = snprintf(command_str,1024,
            "if test -f startscript && test -x startscript;\n"
            "then ./startscript %s %d\n"
            "else ${FINESS}/scripts/startscript %s %d\n"
            "fi", output_dir.c_str(), ndims, output_dir.c_str(), ndims);
    assert_lt(numchars,1023);
    assert_gt(numchars,0);
    int exit_status = system(command_str);
}

