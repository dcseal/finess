#ifndef _INITIALPARAMS_H_
#define _INITIALPARAMS_H_

#include <string>

#include "dog_ini.h"

#include "util.h"


class Params;

extern Params params;

class Params{

public:
    Params(){}
    void init(const std::string& inputFileName);
};



#endif
