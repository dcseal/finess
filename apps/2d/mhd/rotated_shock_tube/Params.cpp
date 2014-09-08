#include "Params.h"

Params params;

void Params::init(const std::string& inputFilename){
    IniDocument ini_doc;
    ini_doc.initFromFile(inputFilename);

}
