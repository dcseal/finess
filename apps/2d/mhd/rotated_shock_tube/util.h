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

template<typename T>
std::string anyToString(const T& a){
    std::stringstream ss;
    ss << a;
    return ss.str();
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

inline std::string& remove_trailing_spaces(std::string& s){
    std::size_t found = s.find_last_not_of(" \t");
    if(found != std::string::npos)
        s.erase(found + 1);
    else
        s.clear();
    return s;
}

inline std::string read_entire_stream(std::istream& is){
    std::string str((std::istreambuf_iterator<char>(is)),
               std::istreambuf_iterator<char>());
    return str;
}

#endif
