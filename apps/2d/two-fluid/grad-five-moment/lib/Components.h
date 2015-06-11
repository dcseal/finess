#define _COMPONENTS_H_

enum ComponentID{
    _rho_i=1,  // rho (ion)
    _M1_i =2,  // 1-momentum (ion)
    _M2_i =3,  // 2-momentum (ion)
    _M3_i =4,  // 3-momentum (ion)
    _N_i  =5,  // energy (ion)
    _rho_e=6,  // rho (elc)
    _M1_e =7,  // 1-momentum (elc)
    _M2_e =8,  // 2-momentum (elc)
    _M3_e =9,  // 3-momentum (elc)
    _N_e  =10, // energy (elc)
    _B1   =11, // 1-magnetic field
    _B2   =12, // 2-magnetic field
    _B3   =13, // 3-magnetic field
    _E1   =14, // 1-electric field
    _E2   =15, // 2-electric field
    _E3   =16, // 3-electric field
    _psi  =17, // B-field cleaning
    _phi  =18, // E-field cleaning
    _entropy_i = 19, // entropy tracking
    _entropy_e = 20, // entropy tracking
};

