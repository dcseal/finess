#include <cstring>        // for strdup
#include "IniDocument.h"  // for ini_doc
#include "OutputDir.h"    // main header file for this class

// -------------------------------------------------------------------------- //
// Single function to read ndims, the number of dimensions from the [dogParams]
// section of parameters.ini.
static int get_ndims()
{
    int ndims;
    const char* s_ndims = ini_doc["dogParams"]["ndims"].c_str();
    if(!sscanf(s_ndims, "%d" ,&ndims))
        eprintf("invalid_value: s_ndims = %s", s_ndims);
    return ndims;
}
// -------------------------------------------------------------------------- //

// convenience accessor
///@brief Name of output directory
///
///- Initialized to 0 so that OutputDir::set_outputdir(...) behave normally;
///- Static member that serves as a global variable.
///  Values can be read by OutputDir::get_outputdir(), which has syncronym ::get_outputdir();
///- OutputDir::parse_arguments sets it to <tt>\<output_some_other_directory\></tt>,
///  if <tt>-o \<output_some_other_directory\> </tt> is provided in command line,
///  and to <tt>output</tt> by default.
char* OutputDir::outputdir   = 0;
///@brief Syncronym to OutputDir::get_outputdir().
const char* get_outputdir() { return OutputDir::get_outputdir(); }

///@brief Sets OutputDir::outputdir to <tt>arg</tt>.
///
///- Implementation detail: would fail if Output::outputdir was not initialized to be 0.
void OutputDir::set_outputdir(const char* arg)
{
    free(outputdir);
    outputdir = strdup(arg);
}

///@brief Default destructor.
OutputDir::~OutputDir()
{ }

///@brief Parses arguments from command line
///
///- Supports <tt>-o \<output_some_other_directory\></tt> and <tt>-d \<debug_level\></tt>;
// -------------------------------------------------------------------------- //
///
///- (Old comments in code, by DS) This function parses user supplied arguments when running the code.
///   e.g., this function parses the section "-o output_some_other_directory"
///   when the user calls "./dog.exe -o output_some_other_directory"
/// 
/// TODO - apparently DEBUG options can be specified from the command line.
/// Because there's no documentation for this, I don't know what all options are
/// available (-DS).
// -------------------------------------------------------------------------- //
///
///- It invokes DebugLevel::set(int) upon parsing <tt><tt>-d \<debug_level\></tt></tt>
///@todo Clean up the code.
void OutputDir::parse_arguments(int argc,char**argv)
{

    // flag used for telling if arguments were sucessfully parsed.
    int show_usage=0;

    // default output directory:
    char* outputdir=(char*)"output";

    // parse arguments
    if(argc>1)
    {
        for(int arg_idx=1; arg_idx<argc; arg_idx++)
        {
            if(argv[arg_idx][0]=='-')
            {

                // "-o" specifies which output directory to use.
                if(argv[arg_idx][1]=='o' && ++arg_idx < argc)
                { outputdir=argv[arg_idx]; }

                // "-d" specifies which debug level to use.   What are the legal
                // options that can be used here? (-DS)
                else if(argv[arg_idx][1]=='d' && ++arg_idx < argc)
                {

                    int debug_level;
                    int success = sscanf(argv[arg_idx],"%d",&debug_level);

                    if(success!=1)
                    { goto show_usage; }

                    DebugLevel::set(debug_level);
                    dprintf("set debug level to %d",DebugLevel::get());
                }
                else{ goto show_usage; }
            }
            else { goto show_usage; }
        }
    }
    set_outputdir(outputdir);
    return;

    show_usage:
    {
        // TODO - add in the debug options as part of this helper statement!
        printf(
                "usage: %s [-o outputdir]\n"
                "  -o outputdir : use outputdir as output directory [default: output]\n",
                argv[0]);
        exit(1);
    }
}
// -------------------------------------------------------------------------- //
