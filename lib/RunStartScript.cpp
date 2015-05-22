#include <string>
#include <stdlib.h> // for system()
#include "assert.h"
#include "debug.h"

#include "IniParams.h"

// Function used to call startscript from $(FINESS)/scripts;
//
// That shell script creates the output directory and copies a few files in
// the output folder
void RunStartScript(std::string parameters_ini_filename, std::string finess_exe_path)
{
    char command_str[1024];
//    const char* get_outputdir();
    std::string output_dir = global_ini_params.get_output_dir();

    // run startscript
    // to create output directory if it does not exist
    // and copy data files to output directory
    int numchars = snprintf(command_str,1024,
            "if test -f startscript && test -x startscript;\n"
            "then ./startscript %s %s %s\n"
            "else ${FINESS}/scripts/startscript %s %s %s\n"
            "fi", output_dir.c_str(),
            parameters_ini_filename.c_str(),
            finess_exe_path.c_str(),
            output_dir.c_str(),
            parameters_ini_filename.c_str(),
            finess_exe_path.c_str());
    assert_lt(numchars,1023);
    assert_gt(numchars,0);
    int exit_status = system(command_str);
}

