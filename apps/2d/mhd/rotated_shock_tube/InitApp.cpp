#include "IniDocument.h"
#include "InitialParams.h"
#include "MHDParams.h"

InitialParams initialParams;
MHDParams mhdParams;

void InitApp(IniDocument& ini_doc)
{
  initialParams.init(ini_doc);
  mhdParams.init(ini_doc);
}
