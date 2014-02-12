#include <stdlib.h> // for system()
#include "assert.h"
#include "debug.h"

///@brief Function used to call @link scripts/startscript @endlink;
///
///- (Old comments in code) That shell script creates the output directory and copies a few files in
/// the output folder
///- It runs the following command on shell
///@verbatim
///if test -f startscript && test -x startscript;
///then ./startscript <output_dir> <ndims>
///else ${FINESS}/scripts/startscript <output_dir> <ndims>
///fi
///@endverbatim
///where <tt>\<output_dir\></tt>=get_outputdir(), and <tt>\<ndims\></tt>=<tt>ndims</tt>
///See also @link scripts/startscript @endlink
///@todo Document @link scripts/startscript @endlink
void RunStartScript(int ndims )
{
    char command_str[1024];
    const char* get_outputdir();

    // run startscript
    // to create output directory if it does not exist
    // and copy data files to output directory
    int numchars = snprintf(command_str,1024,
            "if test -f startscript && test -x startscript;\n"
            "then ./startscript %s %d\n"
            "else ${FINESS}/scripts/startscript %s %d\n"
            "fi", get_outputdir(), ndims, get_outputdir(), ndims);
    assert_lt(numchars,1023);
    assert_gt(numchars,0);
    int exit_status = system(command_str);
}

///@brief Function used to call @link scripts/startscript @endlink.
///
///- See #RunStartScript(int) for what shell commands exactly these two functions execute.
///  Only difference in the current function is <tt>\<output_dir\></tt>=<tt>outputdir</tt>
void RunStartScript(int ndims, const char* outputdir)
{
    char command_str[1024];

    // run startscript
    // to create output directory if it does not exist
    // and copy data files to output directory
    int numchars = snprintf(command_str,1024,
            "if test -f startscript && test -x startscript;\n"
            "then ./startscript %s %d\n"
            "else ${FINESS}/scripts/startscript %s %d\n"
            "fi", outputdir, ndims, outputdir, ndims);
    assert_lt(numchars,1023);
    assert_gt(numchars,0);
    int exit_status = system(command_str);
}
