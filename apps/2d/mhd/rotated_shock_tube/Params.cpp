/// @file Params.cpp
#include "Params.h"

Params params;


template<>
Params::ReconstructionMethod::enum_type stringToAny<Params::ReconstructionMethod::enum_type>(const std::string & s){
    return 
	s == "A" ? Params::ReconstructionMethod::A :
	s == "B" ? Params::ReconstructionMethod::B :
	s == "C" ? Params::ReconstructionMethod::C :
	Params::ReconstructionMethod::DEFAULT;
}


void Params::init(const std::string& inputFilename){
    using std::string;
    IniDocument ini_doc;
    ini_doc.initFromFile(inputFilename);


    string reconstruction_method_str = ini_doc["reconstruction"]["method"];
    if(reconstruction_method_str == "")
    terminate("reconstruction.method is missing.");

    this->reconstruction_method = stringToAny<Params::ReconstructionMethod::enum_type>(reconstruction_method_str);

    if(this->reconstruction_method == Params::ReconstructionMethod::DEFAULT)
        terminate("reconstruction.method must be one of the following: A, B, C.");
    
    string meqn_str = ini_doc["dogParams"]["meqn"];
    if(meqn_str == "")
        terminate("dogParams.meqn is missing.");

    {
        this->meqn = stringToAny<int>(meqn_str);
    }

    if(this->meqn <= 0)
        terminate("dogParams.meqn is expected to be positive.");


}

