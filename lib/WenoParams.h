#ifndef _WENOPARAMS_H_
#define _WENOPARAMS_H_

#include "IniDocument.h"

class WENOParams{
public:
	enum WENOVersion {
		FD, JS, Z
	} weno_version;
	double epsilon;
	int power_param;
	double alpha_scaling;

	void init(IniDocument& ini_doc);
	void append_wenohelp(const char* filename);
};

extern WENOParams wenoParams;

#endif