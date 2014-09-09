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

#include <iostream>
#include <cstdlib>
#include <stdexcept>
//#define TERMINATE_ABORT
#undef TERMINATE_ABORT
inline void terminate(const std::string& message){
#ifdef TERMINATE_ABORT
    std::cerr << message << std::endl;
    std::terminate();
#else
    throw std::runtime_error(message);
#endif
}

#endif
