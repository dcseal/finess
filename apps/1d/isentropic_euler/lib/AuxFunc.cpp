#include <fstream>
#include "tensors.h"
using namespace std;

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(const dTensor1& xpts, 
	     dTensor2& auxvals)
{
  const int numpts=xpts.getsize();
  double x,gamma;
  char buffer[256];
  ifstream read_file("param.data", ios::in);
  
  // Gas consant
  read_file >> gamma;
  read_file.getline(buffer,256);
  read_file.close();
  
  for (int i=1; i<=numpts; i++)
    {
      x = xpts.get(i);
      
      auxvals.set(i,1, gamma );
    }
}
