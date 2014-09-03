#ifndef _INITIALPARAMS_H_
#define _INITIALPARAMS_H_

class IniDocument;
class InitialParams
{  
 public:
  // parameters for initial conditions
  double rhol,unl,utl,u3l,pl,Bnl,Btl,B3l;
  double rhor,unr,utr,u3r,pr,Bnr,Btr,B3r;
  double angle;
  void init(IniDocument& ini_doc);
};
extern InitialParams initialParams;

#endif

