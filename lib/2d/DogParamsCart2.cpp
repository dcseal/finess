#include <stdlib.h> // for exit()
#include <sstream>
#include "DogParamsCart2.h"
#include "IniDocument.h"
#include "debug.h"
using namespace std;

DogParamsCart2 dogParamsCart2;

void DogParamsCart2::set_xlims(double xlow_in, double xhigh_in)
{
    xlow = xlow_in;
    xhigh = xhigh_in;
    dx = (xhigh-xlow)/mx;
}

void DogParamsCart2::set_ylims(double ylow_in, double yhigh_in)
{
    ylow = ylow_in;
    yhigh = yhigh_in;
    dy = (yhigh-ylow)/my;
}

static void not_set_error(string varname)
{
  eprintf("In parameters.ini, section [grid],"
    " you must set the parameter %s\n", varname.c_str());
}

static bool invalid_value(string varname, string val)
{
  eprintf("invalid value: %s = [%s]", varname.c_str(), val.c_str());
  return false;
}

// It does not make sense to define defaults for
// these options, so we do not use the defaults
// mechanism used e.g. in DogParams::init
void DogParamsCart2::init(IniDocument& ini_doc)
{
  if(is_initialized) {
    return;
  }

  IniDocument::Section& grid = ini_doc["grid"];
  //
  string s_mbc   = grid["mbc"  ];
  string s_mx    = grid["mx"   ];
  string s_xlow  = grid["xlow" ];
  string s_xhigh = grid["xhigh"];
  string s_my    = grid["my"   ];
  string s_ylow  = grid["ylow" ];
  string s_yhigh = grid["yhigh"];
  //
  // set defaults (should we read these in from a file?)
  //
  if(s_mbc     .empty()) not_set_error("mbc");
  if(s_mx      .empty()) not_set_error("mx");
  if(s_xlow    .empty()) not_set_error("xlow");
  if(s_xhigh   .empty()) not_set_error("xhigh");
  if(s_my      .empty()) not_set_error("my");
  if(s_ylow    .empty()) not_set_error("ylow");
  if(s_yhigh   .empty()) not_set_error("yhigh");

  // convert strings to numbers
  //
  istringstream is_mbc   (s_mbc   );
  istringstream is_mx    (s_mx    );
  istringstream is_xlow  (s_xlow  );
  istringstream is_xhigh (s_xhigh );
  istringstream is_my    (s_my    );
  istringstream is_ylow  (s_ylow  );
  istringstream is_yhigh (s_yhigh );

  // populate class with parameter data.
  //
  (is_mbc    >> mbc   ) || invalid_value("mbc"   , s_mbc   );
  (is_mx     >> mx    ) || invalid_value("mx"    , s_mx    );
  (is_xlow   >> xlow  ) || invalid_value("xlow"  , s_xlow  );
  (is_xhigh  >> xhigh ) || invalid_value("xhigh" , s_xhigh );
  (is_my     >> my    ) || invalid_value("my"    , s_my    );
  (is_ylow   >> ylow  ) || invalid_value("ylow"  , s_ylow  );
  (is_yhigh  >> yhigh ) || invalid_value("yhigh" , s_yhigh );

  set_xlims(xlow, xhigh);
  set_ylims(ylow, yhigh);

  checkParameters();
  setDerivedParameters();
  reportParameters();

  is_initialized = true;
}

void DogParamsCart2::reportParameters()
{
  printf(  "   === parameters from [grid] ===");
  printf("\n   mx   : %d", mx   );
  printf("\n   my   : %d", my   );
  printf("\n   mbc  : %d", mbc  );
  printf("\n   xlow : %f", xlow );
  printf("\n   xhigh: %f", xhigh);
  printf("\n   ylow : %f", ylow );
  printf("\n   yhigh: %f", yhigh);
  printf("\n   --- parameters derived from [grid] ---");
  printf("\n   dx   : %f ", dx);
  printf("\n   dy   : %f ", dy);
  printf("\n");
}

void DogParamsCart2::setDerivedParameters()
{
    dx = (xhigh-xlow)/mx;
    dy = (yhigh-ylow)/my;
}

void DogParamsCart2::checkParameters()
{
    if(mx <= 0) eprintf("ERROR: mx=%d must be positive.\n", mx);
    if(my <= 0) eprintf("ERROR: my=%d must be positive.\n", my);
    if(mbc < 0) eprintf("ERROR: mbc=%d must be nonnegative.\n", mbc);
    if((xhigh-xlow) <= 0) {
        eprintf("ERROR: xhigh=%d should be greater than"
                " xlow=%d\n", xhigh, xlow);
    }
    if((yhigh-ylow) <= 0) {
        eprintf("ERROR: yhigh=%d should be greater than"
                " ylow=%d\n", yhigh, ylow);
    }

}

// append data to outputdir/qhelp.dat (which is then used by plotting routines)
// TODO - double check that this is the same format used by plotting routines
void DogParamsCart2::append_qhelp(const char* filename)
{
    //dogParams.write_qhelp(filename);
    FILE* file = fopen(filename,"a");
    fprintf(file,"%16d : mx\n", mx);
    fprintf(file,"%16d : my\n", my);
    fprintf(file,"%16.8e : xlow\n", xlow);
    fprintf(file,"%16.8e : xhigh\n",xhigh);
    fprintf(file,"%16.8e : ylow\n", ylow);
    fprintf(file,"%16.8e : yhigh\n",yhigh);
    fclose(file);
}
