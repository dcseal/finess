#include "dogdefs.h"
#include "MHDParams.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// Function that is called after initial conditions are set
void AfterQinit(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q)
{

  const int melems = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   maux = aux.getsize(2);  
  
  // Output parameters to file in outputdir
  /*
  char mhdhelp[200];
  strcpy(mhdhelp,solver.get_outputdir());
  strcat(mhdhelp,"/mhdhelp.dat");
  mhdParams.write_mhdhelp(mhdhelp);
  */
}
