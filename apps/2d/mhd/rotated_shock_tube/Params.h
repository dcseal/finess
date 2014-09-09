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

public:
    struct ReconstructionMethod{
    	enum enum_type {A, B, C, DEFAULT};
    };
private:    
    ReconstructionMethod::enum_type reconstruction_method;
public:
    inline ReconstructionMethod::enum_type get_reconstruction_method(){
        return this->reconstruction_method;
    }
private:
    int meqn;
public:
    inline int get_meqn(){
	return this->meqn;
    }
private:
    
};



#endif
