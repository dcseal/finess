#include "tensors.h"
#include "IniParams.h"

void GridSetup( dTensor3& node, dTensor2& prim_vol)
{

    const double xlow = global_ini_params.get_xlow();
    const double dx   = global_ini_params.get_dx();

    const double ylow = global_ini_params.get_ylow();
    const double dy   = global_ini_params.get_dy();

// TODO - ??? - mnodes = mx*my, or mnodes = (mx,my) ???
//  const int        = node.getsize(1);
//  const int melems = prim_vol.getsize();

//      // Set variable "node"
//  #pragma omp parallel for
//      for (int i=1; i<=mnodes; i++)
//      {  node.set(i,1, xlow+(double(i)-1.0e0)*dx );  }

//      // Set grid cell volumes
//  #pragma omp parallel for
//      for (int j=1; j<=melems; j++)
//      {  prim_vol.set(j, dx);  }

}
