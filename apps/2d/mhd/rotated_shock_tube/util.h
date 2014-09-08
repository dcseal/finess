#include <sstream>
template<typename T>
T stringToAny(const std::string& s){
    T ret;
    std::istringstream iss(s);
    iss >> ret;
    return ret;
}


