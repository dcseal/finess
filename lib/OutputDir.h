#ifndef _OUTPUT_DIR_H__ 
#define _OUTPUT_DIR_H__ 

#include "debug.h"

// -------------------------------------------------------------------------- //
// Used to create output directory and copy a few files (Makefiles,
// parameters.ini, ... ) into the output folder.
// -------------------------------------------------------------------------- //
///@brief A class that parses command line options.
///
///- Poorly named.  Should have been called CommandLineOptions.
///  See #OutputDir::parse_arguments(int, char**) for details.
///- Poorly organized.  #parse_arguments(int, char**) does all the job, 
///  and it does so by setting global variables and returns <tt>void</tt> (!!!)
///@todo Merge this with #IniDocument into a common program-options class.
class OutputDir
{

    // saving and initialization
    public:

        // ------------------------------------------------------------------ //
        // parser to read input arguments.  Current valid options are 
        // "-o outputdir" and "-d X", for an output directory and debugging
        // level.
        // ------------------------------------------------------------------ //
        static void parse_arguments(int argc, char**argv);
        ///@brief Returns OutputDir::outputdir
        static const char* get_outputdir() {return outputdir;}
        virtual ~OutputDir();

    private:
        static char* outputdir;

    private:
        static void set_outputdir(const char*);

};

// convenience global accessor
const char* get_outputdir();

#endif
