#include <sstream>
#include <fstream>
#include <string>
#include "stdlib.h"

#include "IniParams.h"
#include "StateVars.h"

using namespace std;

// This routine is called once per simulation.  The purpose of this routine is
// to provide access to the plotting routines provided from DoGPack.
//
// Make sure that the global IniParams has already been called before calling
// this function.
//
// See: $FINESS/viz/matlab/plotfin1.m for where this gets parsed in the
// plotting routines.
void WriteQhelp( )
{

    string outputdir = global_ini_params.get_output_dir();

    ostringstream fname;
    fname << outputdir << "/" << "qhelp.dat";

    // Open file -- qhelp.dat
    ostringstream fname1;
    fname1 << outputdir << "/" << "qhelp.dat";
    FILE* file = fopen( fname1.str().c_str(), "w" );

    fprintf(file,"%16d   : ndims   \n", global_ini_params.get_ndims() );
    fprintf(file,"%16d   : meqn    \n", global_ini_params.get_meqn()  );
    fprintf(file,"%16d   : maux    \n", global_ini_params.get_maux()  );
    fprintf(file,"%16d   : nout    \n", global_ini_params.get_nout()  );
    fprintf(file,"%16d   : mx      \n", global_ini_params.get_mx()    );
    fprintf(file,"%16.8e : xlow    \n", global_ini_params.get_xlow()  );
    fprintf(file,"%16.8e : xhigh   \n", global_ini_params.get_xhigh() );
    fclose(file);

}
