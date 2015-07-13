#ifndef _MAXWELL_H_
#define _MAXWELL_H_
#include "math.h"
#include "debug.h"

enum MaxwellComponentID{
    M_B1   = 1, // 1-magnetic field
    M_B2   = 2, // 2-magnetic field
    M_B3   = 3, // 3-magnetic field
    M_E1   = 4, // 1-electric field
    M_E2   = 5, // 2-electric field
    M_E3   = 6, // 3-electric field
    M_psi  = 7, // B-field cleaning
    M_phi  = 8, // E-field cleaning
};

enum MaxwellEigID {
    M_m_B2_E3  = 1,
    M_p_B2_E3  = 2,
    //
    M_m_B3_E2  = 3,
    M_p_B3_E2  = 4,
    //
    M_m_B1_psi = 5,
    M_p_B1_psi = 6,
    //
    M_m_E1_phi = 7,
    M_p_E1_phi = 8,
};

enum SymMaxwellComponentID{
    SM_B1   = 1, // 1-magnetic field
    SM_B2   = 2, // 2-magnetic field
    SM_E3   = 3, // 3-electric field
    SM_psi  = 4, // B-field cleaning
};

enum SymMaxwellEigID {
    SM_m_B2_E3  = 1,
    SM_p_B2_E3  = 2,
    //
    SM_m_B1_psi = 3,
    SM_p_B1_psi = 4,
};

/*
 * The routines below here were pulled from DoGPack.  I am not sure what they
 * were used for. -DS 7/9/2015

// This function oscillates between 0 and 1 with
// a period of 4.  The definition is:
//
// Let theta = t*eps.
// f(theta) = f(theta mod 4) and
// f(theta) = 0 if 0 <= theta <= 1,
// f(theta) increases smoothly from 0 to eps as theta increases from 1 to 2
// f(theta) = 1 if 2 <= theta <= 3,
// f(theta) decreases smoothly from 1 to 0 as theta increases from 3 to 4
//
inline double modulate_source(double t, double eps)
{
  double theta = eps*t;
  // return eps*pow(sin(M_PI_2*theta),2);
  double theta_mod_4 = fmod(theta,4);
  switch(int(theta_mod_4))
  {
    case 0: return 0.;
    case 1:
    {
      const double s = sin(M_PI_2*(theta_mod_4-1.));
      return eps*s*s;
    }
    case 2: return eps;
    case 3:
    {
      const double s = sin(M_PI_2*(theta_mod_4-2.));
      return eps*s*s;
    }
    default:
      invalid_value_error(theta);
  }
}

inline double eps_Bcp_func(double psi, double time)
{
  return psi*modulate_source(time,maxwellParams.get_eps_Bcp());
}

inline double eps_Ecp_func(double phi, double time)
{
  return phi*modulate_source(time,maxwellParams.get_eps_Ecp());
}
*/

#endif
