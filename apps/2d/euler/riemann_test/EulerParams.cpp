#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "EulerParams.h"
#include "debug.h"

void EulerParams::init(IniDocument& ini_doc)
{
  string section_label = "eulerParams";
  IniDocument::Section& ini_sec = ini_doc[section_label];

  //
  // Read-in information into strings from file
  // "parameters.ini" from section [mhdparams]
  //
  string s_gamma             = ini_sec["gamma"         ];
  string s_OPT               = ini_sec["OPT"           ];
  
  //
  // Convert read-in strings to numerical values
  //
  istringstream is_gamma           (s_gamma         );
  istringstream is_OPT             (s_OPT           );

  //
  // Store information in global parameter values
  //
  gamma = 1.6666666666666667;
  OPT = 1;

  is_gamma >> gamma;  
  is_OPT >> OPT;

  //
  // Output values to screen
  //
  cout << "   === parameters from [" << section_label << "] ===" << endl;
  cout << "   gamma             :  " << gamma           << endl;
  cout << "   OPT               :  " << OPT           << endl;
  cout << endl;
}

void EulerParams::write_eulerhelp(const char* filename)
{
  FILE* file = fopen(filename,"w");
  fprintf(file,"%24.16e\n", gamma);
  fprintf(file,"%d\n", OPT);
  fclose(file);
}
