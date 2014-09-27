#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "IniParams.h"
#include "IniParams.h"
using namespace std;

// Template used for outputting extra information.
// Applicaitons can link out this template if desired.
void Output_Extra(
            const dTensorBC2& aux,
            const dTensorBC2& q,
            double t,
            int nframe )
{

    // See $(FINESS)/lib/1d/Output.cpp to see where this gets called.
    return;

}
