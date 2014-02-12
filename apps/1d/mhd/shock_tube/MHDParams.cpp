#include<assert.h>
#include<iostream>
#include<sstream>
#include<stdlib.h> // for exit()
#include <cmath>
#include "dog_ini.h"
#include "MHDParams.h"

void MHDParams::init(IniDocument& ini_doc)
{
    string section_label = "mhdparams";
    IniDocument::Section& ini_sec = ini_doc[section_label];

    //
    // Read-in information into strings from file
    // "parameters.ini" from section [mhdparams]
    //
    string s_gamma             = ini_sec["gamma"         ];

    //
    // Convert read-in strings to numerical values
    //
    istringstream is_gamma           (s_gamma         );

    //
    // Store information in global parameter values
    //
    gamma = 1.6666666666666667;

    is_gamma >> gamma;  

    //
    // Output values to screen
    //
    cout << "   === parameters from [" << section_label << "] ===" << endl;
    cout << "   gamma             :  " << gamma           << endl;
    cout << endl;
}

void MHDParams::write_mhdhelp(const char* filename)
{
    FILE* file = fopen(filename,"w");
    fprintf(file,"%24.16e\n", gamma);
    fclose(file);
}
