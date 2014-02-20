///@file apps/1d/euler/harten_shock_tube/main.cpp
#include "dogdefs.h"

// =========================================================================
//
//  Copyright J.A. Rossmanith
//
//  This software is made available for research and instructional use only.
//  You may copy and use this software without charge for these non-commercial
//  purposes, provided that the copyright notice and associated text is
//  reproduced on all copies.  For all other uses (including distribution of
//  modified versions), please contact the author at the address given below.
//
//  *** This software is made available "as is" without any assurance that it
//  *** will work for your purposes.  The software may in fact have defects, so
//  *** use the software at your own risk.
//
//  -------------------------------
//  DoGPack
//  -------------------------------
//
//    Lead Developer:  
//             James Rossmanith
//             Iowa State University
//             Department of Mathematics
//             396 Carver Hall
//             Ames, IA 50011
//             rossmani@iastate.edu
// =========================================================================

///@brief Calls #main_global(int, char*[]) and returns its return value.
///
//
int main(int argc, char* argv[])
{
 
  int m;
  int main_global(int argc, char* argv[]);
  m = main_global(argc,argv);
  
  return m;
}
