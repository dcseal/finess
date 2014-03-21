#ifndef _WENOPARAMS_H_
#define _WENOPARAMS_H_

// TODO - write some documention for this module.  For example, how does a
// user access this stuff?  See also statements?

#include "IniDocument.h"

///@brief Class that stores WENO Parameters, which are read from parameters.ini file.
///
///Usage: see #wenoParams
class WENOParams{
public:
	///@brief Reconstruction method.
	///
	///- <tt>FD</tt>  Finite difference (??)
	///- <tt>JS</tt>  Jiang and Shu
	///- <tt>Z</tt>   Z
	enum WENOVersion {
		FD, JS, Z
	} weno_version;
	///@brief  @f$ \epsilon @f$ in equation (2.16), JS1996
	double epsilon;
	///@brief  @f$ p @f$ in equation (2.16), JS1996
	double power_param;
	///@brief  What is this??  Why here??
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
