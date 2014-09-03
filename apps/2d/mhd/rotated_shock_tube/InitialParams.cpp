#include<assert.h>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "InitialParams.h"

void InitialParams::init(IniDocument& ini_doc)
{
  string section_label = "initialparams";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  //
  // Read-in information into strings from file
  // "parameters.ini" from section [initialparams]
  //
  string s_rhol  = ini_sec["rhol" ];
  string s_unl   = ini_sec["unl"  ];
  string s_utl   = ini_sec["utl"  ];
  string s_u3l   = ini_sec["u3l"  ];
  string s_pl    = ini_sec["pl"   ];
  string s_Bnl   = ini_sec["Bnl"  ];
  string s_Btl   = ini_sec["Btl"  ];
  string s_B3l   = ini_sec["B3l"  ];

  string s_rhor  = ini_sec["rhor" ];
  string s_unr   = ini_sec["unr"  ];
  string s_utr   = ini_sec["utr"  ];
  string s_u3r   = ini_sec["u3r"  ];
  string s_pr    = ini_sec["pr"   ];
  string s_Bnr   = ini_sec["Bnr"  ];
  string s_Btr   = ini_sec["Btr"  ];
  string s_B3r   = ini_sec["B3r"  ];

  string s_angle = ini_sec["angle"];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream is_rhol  (s_rhol );
  istringstream is_unl   (s_unl  );
  istringstream is_utl   (s_utl  );
  istringstream is_u3l   (s_u3l  );
  istringstream is_pl    (s_pl   );
  istringstream is_Bnl   (s_Bnl  );
  istringstream is_Btl   (s_Btl  );
  istringstream is_B3l   (s_B3l  );

  istringstream is_rhor  (s_rhor );
  istringstream is_unr   (s_unr  );
  istringstream is_utr   (s_utr  );
  istringstream is_u3r   (s_u3r  );
  istringstream is_pr    (s_pr   );
  istringstream is_Bnr   (s_Bnr  );
  istringstream is_Btr   (s_Btr  );
  istringstream is_B3r   (s_B3r  );

  istringstream is_angle (s_angle);

  //
  // Store information in global parameter values
  //
  is_rhol >> rhol;
  is_unl  >> unl;
  is_utl  >> utl;
  is_u3l  >> u3l;
  is_pl   >> pl;
  is_Bnl  >> Bnl;
  is_Btl  >> Btl;
  is_B3l  >> B3l;
  
  is_rhor >> rhor;
  is_unr  >> unr;
  is_utr  >> utr;
  is_u3r  >> u3r;
  is_pr   >> pr;
  is_Bnr  >> Bnr;
  is_Btr  >> Btr;
  is_B3r  >> B3r;

  is_angle >> angle;
  //
  // Output values to screen
  //
  cout << scientific;
  cout << setprecision(6);
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   rhol   :  " << setw(14) << rhol << endl;
  cout << "   unl    :  " << setw(14) << unl << endl;
  cout << "   utl    :  " << setw(14) << utl << endl;
  cout << "   u3l    :  " << setw(14) << u3l << endl;
  cout << "   pl     :  " << setw(14) << pl << endl;
  cout << "   Bnl    :  " << setw(14) << Bnl << endl;
  cout << "   Btl    :  " << setw(14) << Btl << endl;
  cout << "   B3l    :  " << setw(14) << B3l << endl;
  cout << endl;
  cout << "   rhor   :  " << setw(14) << rhor << endl;
  cout << "   unr    :  " << setw(14) << unr << endl;
  cout << "   utr    :  " << setw(14) << utr << endl;
  cout << "   u3r    :  " << setw(14) << u3r << endl;
  cout << "   pr     :  " << setw(14) << pr << endl;
  cout << "   Bnr    :  " << setw(14) << Bnr << endl;
  cout << "   Btr    :  " << setw(14) << Btr << endl;
  cout << "   B3r    :  " << setw(14) << B3r << endl;
  cout << endl;
}
