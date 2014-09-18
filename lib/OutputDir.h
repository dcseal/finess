#error "Stop using OutputDir!"
#ifndef _OUTPUT_DIR_H__ 
#define _OUTPUT_DIR_H__ 

#include "debug.h"

// -------------------------------------------------------------------------- //
// Used to create output directory and copy a few files (Makefiles,
// parameters.ini, ... ) into the output folder.
// -------------------------------------------------------------------------- //
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
