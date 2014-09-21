#include "dogdefs.h"
#include <string.h>          // For strcpy and strcat (some compilers don't need this)
#include "IniParams.h"

const char* get_outputdir();

// Function that is called after initial condition
void AfterQinit(dTensorBC3& aux, dTensorBC3& q)
{

    // Output parameters to file in outputdir
    char eulerhelp[200];
    strcpy( eulerhelp, get_outputdir() );
    strcat( eulerhelp, "/eulerhelp.dat");
    eulerParams.write_eulerhelp( eulerhelp );

}
