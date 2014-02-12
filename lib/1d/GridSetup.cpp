#include "DogParamsCart1.h"
#include "tensors.h"

///@brief Sets up arrays containing cell edges and volume of each element.
///
///- Currently, nodes (i.e. cell edges) are uniform, starting from <xlow>, with step <dx>
///  where <xlow> and <dx> are extracted from #ini_doc by methods in class #dogParamsCart1.
///@param node Is set to contain cell edges
///@param prim_vol Is set to contain the volume (length) of each element.
void GridSetup( dTensor2& node, dTensor1& prim_vol)
{
    const double xlow = dogParamsCart1.get_xlow();
    const double dx   = dogParamsCart1.get_dx();

    const int mnodes = node.getsize(1);
    const int melems = prim_vol.getsize();

    // Set variable "node"
#pragma omp parallel for
    for (int i=1; i<=mnodes; i++)
    {  node.set(i,1, xlow+(double(i)-1.0e0)*dx );  }

    // Set grid cell volumes
#pragma omp parallel for
    for (int j=1; j<=melems; j++)
    {  prim_vol.set(j, dx);  }

}
