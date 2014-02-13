#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart1.h"
using namespace std;

// Template used for outputting extra information.
// Applicaitons can link out this template if desired.
///@brief Outputs extra information.  See #Output(...) for format information.
///
///- Currently, does nothing.
void Output_Extra(const dTensor2& node, 
            const dTensorBC2& aux,
            const dTensorBC2& q,
            double t,
            int nframe,
            string outputdir)
{

    // See $(FINESS)/lib/1d/Output.cpp to see where this gets called.
    return;

}
