#error "Stop using WenoParams!"
#ifndef _WENOPARAMS_H_
#define _WENOPARAMS_H_

#include "IniDocument.h"

///@brief Class that stores WENO Parameters, which are read from parameters.ini file.
///
///Usage: see #wenoParams
class WENOParams{
public:
    ///@brief Reconstruction method.
    ///
    ///- <tt>FD</tt>  Conservative reconstruction (linear weights)
    ///- <tt>JS</tt>  Jiang and Shu smoothness indicators
    ///- <tt>Z</tt>   Z (WENO-Z - improved order of accuracy)
    enum WENOVersion {
        FD, JS, Z
    } weno_version;

    ///@brief  @f$ \epsilon @f$ in equation (2.16), JS1996
    double epsilon;

    ///@brief  @f$ p @f$ in equation (2.16), JS1996
    double power_param;

    ///@brief  @f$ \alpha @f$ scaling used in the Lax-Friedrich's flux splitting
    double alpha_scaling;

    ///@brief Populate the class member variables with data from section [weno] of <tt>ini_doc</tt>
    ///
    void init(IniDocument& ini_doc);
    ///@brief ???
    //All the other xxxParams has write_xxxhelp member function, so just copied one here
    void append_wenohelp(const char* filename);
};

extern WENOParams wenoParams;

#endif
