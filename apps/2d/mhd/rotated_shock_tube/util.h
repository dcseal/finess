#ifndef _UTIL_H_
#define _UTIL_H_

#include <sstream>
template<typename T>
T stringToAny(const std::string& s){
    T ret;
    std::istringstream iss(s);
    iss >> ret;
    return ret;
}

#include <cstdio>
inline bool existFile(const std::string& filename){
    using std::FILE;
    using std::fopen;
    using std::fclose;
    if(FILE *file = fopen(filename.c_str(), "r")){
	fclose(file);
	return true;
    } else{
	return false;
    }
}


#endif
