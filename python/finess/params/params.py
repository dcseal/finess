"""Provides framework for generating IniParams, a C++ class for
interpreting .ini config file, and global_ini_params, a C++ global
variable of type IniParams.

See files under _templates/ for sample usage of this framework.

Example:
If your app is 2D, and needs only those parameters common to all 2D
apps, then you need only
- copy _templates/dim2/simple/generate_iniparams.py to
you app directory, and 
- add #include "IniParams.h" to all the source files you want to use
  the parameters.

"""


# TODO Fully document this module
# For now, 
#   ...
#   [section]
#   name = value
#   ...
#   is first read into C++ local variable variable_name_str,
#      then type-converted via stringToAny<T> template function,
#      and stored into global_ini_params.variable_name.
#
#   stringToAny uses stringstream, and must undergo a
#   template-specialization in order to interpret non-primitive types.
#
#   Currently, the only non-primitive type is EnumParameterType's.
#   In the argument to its constructor, string_to_enumerator_dict
#   refers to the map from what you want to see in .ini, to what you
#   want to see in C++ source, e.g. {"Runge-Kutta": "RK",
#   "Lax-Wendroff": "LxW"}
#
#   A trick in C++ code is used here to isolate namespace of the
#   enumerators (i.e. the namespace of RK, LxW in above example).
#   This namespace is called enum_scope_name in the Python code.
#   For example, if enum_scope_name is "TimeSteppingMethod", then full
#   name of RK above is IniParams::TimeSteppingMethod::RK, which is of
#   type IniParams::TimeSteppingMethod::enum_type.  This is to
#   circumvent a limitation in C++98.  Otherwise a C++11 'scoped enum'
#   could be used.
#   
#   Accessor was meant to provide methods of getting an attribute in a
#   different manner from global_ini_params.get_xxx().  However, since
#   all C++ code can be modified to use this standard get_xxx method,
#   Accessor seems redundant now.
#
#   There are two places where a user needs to insert C++ code
#   directly as the argument to a Python call:
#
#   - Defining expression of a DerivedParameter
#   - Non-canonical Check.
#
#   Usually, a user only needs to know how to use the following:
#   - Constructors of various Parameter classes
#   - Check classes
#   - append_pac_from_module
#   - write_to_header_cpp
#
#   The usage of items in this list should be clear from the templates
#   in _templates, and from the source code of finess.viz.dim2



from __future__ import absolute_import

class_name = "IniParams"
global_variable_name = "global_ini_params"




def terminate_on_missing(variable_name, section, name):
    return """
    if(%(variable_name)s_str == "")
        terminate("%(section)s.%(name)s is missing.");
    """ % \
    {"variable_name" : variable_name,
     "section": section, 
     "name": name}



def generate_default_on_missing(default_value):
    def default_on_missing(variable_name, section, name):
        return """
    if(%(variable_name)s_str == ""){
        %(variable_name)s_str = "%(default_value)s";
	this->ini_doc["%(section)s"]["%(name)s"] = "%(default_value)s";
    }
        """ % \
        {"variable_name" : variable_name,
         "section": section, 
         "name": name,
         "default_value": \
	     str(default_value) if type(default_value) != bool \
	         else str(default_value).lower()}
    return default_on_missing



class ParameterType:
    def __init__(self, type_string):
        if type(type_string) != str:
            raise ValueError("type_string should be a string.")                    
                             
        canonical_types = {"int", "double", "std::string", "bool"}
        if not type_string in canonical_types:
            raise ValueError("type_string should be one of " + str(canonical_types))
        self.type_string = type_string
        
    def generate_declaration(self, variable_name):
        declaration = """private:
        %(type_string)s %(variable_name)s;
        """ % \
        {"type_string": self.type_string,
         "variable_name": variable_name}
        
        return declaration
    
    def generate_stringTo_code(self, variable_name):
        return """
    this->%(variable_name)s = stringToAny<%(type_string)s>(%(variable_name)s_str);
        """ % \
        {"variable_name": variable_name,
         "type_string": self.type_string
         }
    
    def generate_parsing(self, variable_name, section, name, 
                         missing_handler = terminate_on_missing):
        parsing = """
    string %(variable_name)s_str = this->ini_doc["%(section)s"]["%(name)s"];
    %(missing_handling)s
    %(stringTo_code)s
    """ % \
        {"variable_name": variable_name,
         "section": section,
         "name": name,
         "missing_handling": missing_handler(variable_name, section, name),
         "stringTo_code": self.generate_stringTo_code(variable_name)}
        return parsing
    



class EnumParameterType(ParameterType):
    def __init__(self, enum_scope_name, string_enumerator_dict):
        if type(enum_scope_name) != str:
            raise ValueError("enum_scope_name should be a string.")
        if type(string_enumerator_dict) != dict:
            raise ValueError("string_enumerator_dict should be a dict.")
        
        self.enum_scope_name = enum_scope_name

        self.string_enumerator_dict = string_enumerator_dict

    @property
    def full_enum_scope_name(self):
        return class_name + "::" + self.enum_scope_name
    
    @property
    def type_string(self):
        return """%(class_name)s::%(enum_scope_name)s::enum_type"""% \
               {"class_name": class_name,
	        "enum_scope_name": self.enum_scope_name}

        

    def generate_type_definition(self):
        return """public:
        struct %(enum_scope_name)s{
            enum enum_type {%(enumerator_list)s, DEFAULT};
        };
        """ % \
        {"enum_scope_name": self.enum_scope_name,
         "enumerator_list": ", ".join(self.string_enumerator_dict.values())}
    
    def generate_stringToAny_specialization(self):
        return """template<>
inline %(type_string)s stringToAny<%(type_string)s>(const std::string& s){
    return
%(case_code)s
        
}
        """ % \
        {"type_string": self.type_string,
         "full_enum_scope_name": self.full_enum_scope_name,
         "case_code":
             "\n".join(["""        s == "%(string)s" ? %(full_enum_scope_name)s::%(enumerator)s :""" % \
                        {"string": item[0],
                         "enumerator": item[1],
                         "full_enum_scope_name": self.full_enum_scope_name} \
                        for item in self.string_enumerator_dict.items()]    + 
                       ["        " + self.full_enum_scope_name + "::DEFAULT;"])
         }
    
    def generate_parsing(self, variable_name, section, name, 
                     missing_handler = terminate_on_missing):
        return ParameterType.generate_parsing(self, variable_name, section, name, missing_handler) +                """
    if(this->%(variable_name)s == %(full_enum_scope_name)s::DEFAULT)
        terminate("%(section)s.%(name)s should be one of the following: %(string_list)s.");
               """ % \
               {"variable_name": variable_name,
                "full_enum_scope_name": self.full_enum_scope_name,
                "section": section,
                "name": name,
                "string_list": ", ".join(self.string_enumerator_dict.keys())}


class Parameter:
    def __init__(self, variable_name, section, name, type_,
                 default_value = None):
        assert type(variable_name) == str
        assert type(section) == str
        assert type(name) == str
        self.variable_name = variable_name
        self.section = section
        self.name = name
        
        if type(type_) == str:
            self.type_= ParameterType(type_)
        elif type(type_) == type(self):
            self.type_ = type_
        else:
            raise ValueError("type_ should be a string, or a ParameterType object.")
	
	self.missing_handler = terminate_on_missing if default_value == None else \
	                       generate_default_on_missing(default_value)

    
    def get_declaration_code(self):
        return self.type_.generate_declaration(self.variable_name)
    
    def get_defining_code(self):
        return self.type_.generate_parsing(self.variable_name, 
	                                   self.section, self.name,
	                                   missing_handler = \
					       self.missing_handler)


class DerivedParameter(Parameter):
    def __init__(self, variable_name, type_, defining_expression_in_cpp):
        assert type(variable_name) == str
               
        if type(type_) == str:
            self.type_= ParameterType(type_)
        elif type(type_) == type(self):
            self.type_ = type_
        else:
            raise ValueError("type_ should be a string, or a ParameterType object.")
        
        assert type(defining_expression_in_cpp) == str
        
        self.variable_name = variable_name
        self.defining_expression_in_cpp = defining_expression_in_cpp
        
    def get_defining_code(self):
        return "    this->%(variable_name)s = %(defining_expression)s;"  %                {"variable_name": self.variable_name,
                "defining_expression": self.defining_expression_in_cpp}

        
class Accessor:
    def __init__(self, parameter, access_by_name = None):
        assert isinstance(parameter, Parameter)
        
        self.parameter = parameter
        
        assert access_by_name == None or type(access_by_name) == str
        self.access_by_name = parameter.variable_name if access_by_name == None else access_by_name
        
    def get_accessor_code(self):
        accessor = """public:
        inline %(type_name)s get_%(access_by_name)s(){
            return this->%(variable_name)s;
        }
        """ % \
        {"type_name": self.parameter.type_.type_string,
         "access_by_name": self.access_by_name,
         "variable_name": self.parameter.variable_name}
        
        return accessor
        

class Check:
    def __init__(self, cpp_code = ""):
        self.cpp_code = cpp_code
    
    def get_cpp_code(self):
        return self.cpp_code


def generate_compare_predicate_Check(compare_op):
    assert compare_op in ["<", "<=", ">", ">=", "=="]
    class CompareCheck(Check):
        def __init__(self, parameter, value):
            assert isinstance(parameter, Parameter)
            self.parameter = parameter
            self.value = value
            self.compare_op = compare_op
        def get_cpp_code(self):
            return """
    if(!(this->%(variable_name)s %(compare_op)s %(value)s))
        terminate("%(variable_name)s should be %(compare_op)s %(value)s.");
""" % \
               {"variable_name": self.parameter.variable_name,
                "compare_op": compare_op,
                "value": self.value}                
            
    CompareCheck.__doc__ = """Check class for "%(compare_op)s".""" % {"compare_op": compare_op}
    return CompareCheck
        

CheckGreaterThan = generate_compare_predicate_Check(">")
CheckLessThan = generate_compare_predicate_Check("<")
CheckGreaterEqual = generate_compare_predicate_Check(">=")
CheckLessEqual = generate_compare_predicate_Check("<=")


class CheckOneOf(Check):
    def __init__(self, parameter, list_):
        assert isinstance(parameter, Parameter)
        assert parameter.type_.type_string in ["int"]
        assert type(list_) in [list, set]
        assert all(map(lambda x: type(x) == int, list_))
        
        self.parameter = parameter
        self.list_ = list_
    
    def get_cpp_code(self):
        return "    if(" +                " &&\n       ".join(["%(variable_name)s != %(item)d" %                                    {"variable_name": self.parameter.variable_name,
                                    "item": item} 
                                   for item in self.list_])  + \
               """)
            terminate("%(variable_name)s is not one of %(list_string)s.");
""" % {"variable_name": self.parameter.variable_name, "list_string": str(self.list_)}


def generate_header_cpp(parameters, accessors, checks):
    header_filename = class_name + ".h"
    cpp_filename = class_name + ".cpp"

    import inspect, os
    current_python_filename = inspect.getfile(inspect.currentframe())
    header_comments_file = """/// @file %(header_filename)s
/// Generated by %(current_python_filename)s
    """ % \
    {"header_filename": header_filename, "current_python_filename": current_python_filename}
    
    cpp_comments_file = """/// @file %(cpp_filename)s
/// Generated by %(current_python_filename)s
""" % \
{"cpp_filename": cpp_filename, "current_python_filename": current_python_filename}

    header_guard_macro_name = "_" + class_name.upper() + "_H_"
    
    header_comments_usage = """
// This file, %(header_filename)s, and %(cpp_filename)s, together define a class %(class_name)s,
// and a global variable %(global_variable_name)s.
//
// To use %(global_variable_name)s,
// 1. Add #include "%(header_filename)s" to all the files that use %(global_variable_name)s;
// 2. Add %(global_variable_name)s.init(<parameter_filename>); to main function, where
//    <parameter_filename> is to be replaced by the actual parameter filename;
// 3. Make sure -I. is in CXXFLAGS of Makefile;
// 4. Make sure %(cpp_filename)s is compiled and linked.
""" % \
{"header_filename": header_filename,
 "cpp_filename": cpp_filename,
 "class_name": class_name,
 "global_variable_name": global_variable_name}
    
    header_head = """#ifndef %(header_guard_macro_name)s
#define %(header_guard_macro_name)s
%(header_comments_file)s

#include <cmath>
#include <string>
#include "util.h"
#include "IniParser.h"

class %(class_name)s;
extern %(class_name)s %(global_variable_name)s;

""" % \
{"header_comments_file": header_comments_file,
 "header_guard_macro_name": header_guard_macro_name,
 "class_name": class_name,
 "global_variable_name": global_variable_name }

    header_tail = """
#endif
"""
    for p in parameters:
        if isinstance(p.type_, EnumParameterType):
	    p.type_.class_name = class_name

    types_requiring_definitions =         [t for t in set([p.type_ for p in parameters])
             if getattr(t, "generate_type_definition", None) != None]
   
    class_declaration_contents =         "// Type definitions \n" +         "\n".join([t.generate_type_definition() for t in types_requiring_definitions]) +         "\n\n// Member variables declarations\n" +         "\n".join([p.get_declaration_code() for p in parameters]) +         "\n\n// Accessor definitions\n" +         "\n".join([a.get_accessor_code() for a in accessors])
        
        
    types_requiring_stringToAny_specialization =         [t for t in set([p.type_ for p in parameters])
             if getattr(t, "generate_stringToAny_specialization", None) != None]
    stringToAny_specializations_code =         "//stringToAny specializations\n" +         "\n".join([t.generate_stringToAny_specialization() for t in types_requiring_stringToAny_specialization])
    
    class_declaration = """
class %(class_name)s{
public:
    void init(const std::string& inputFilename);
private:
    IniParser::ini_doc_type ini_doc;
public:
    std::string ini_doc_as_string(){
        return IniParser::ini_doc_as_string(this->ini_doc);
    }
%(class_declaration_contents)s
};
""" % \
{"class_name": class_name, "class_declaration_contents": class_declaration_contents}

    cpp_head = """%(cpp_comments_file)s
#include "%(header_filename)s"

#include "util.h"
#include "IniParser.h"

#include <fstream>
#include <string>

%(class_name)s %(global_variable_name)s;

""" % \
{"cpp_comments_file": cpp_comments_file,
 "header_filename": header_filename,
 "class_name": class_name,
 "global_variable_name": global_variable_name}

    cpp_init_method_contents = "// Defining code for member variables\n\n" +                                "\n".join([p.get_defining_code() for p in parameters]) +                                "\n// Checks\n\n" +                                "\n".join([c.get_cpp_code() for c in checks])
            
    cpp_init_method = """
void %(class_name)s::init(const std::string& inputFilename){
    using std::string;
    
    IniParser parser;
    {
        std::ifstream ifs(inputFilename.c_str());
	std::string ini_file_content = read_entire_stream(ifs);
	int parse_return_value = parser.parse(ini_file_content);
	if(parse_return_value != 0)
	    terminate("Error parsing " + inputFilename + ": line #" +
	              anyToString(parse_return_value));
    }

    this->ini_doc = parser.get_ini_doc();

%(cpp_init_method_contents)s

}

""" % \
{"class_name":class_name, "cpp_init_method_contents": cpp_init_method_contents}

    header_code = header_head + class_declaration + stringToAny_specializations_code + header_tail
    cpp_code = cpp_head + cpp_init_method
    
    return header_filename, header_code, cpp_filename, cpp_code
    


def append_pac_from_module(pac, module):
    parameters, accessors, checks = pac

    parameters += module.parameter_list
    accessors += module.accessor_list
    checks += module.check_list
    

def write_to_header_cpp(parameter_list, accessor_list, check_list):
    from finess.params import generate_header_cpp
    header_filename, header_code, cpp_filename, cpp_code = generate_header_cpp(parameter_list, accessor_list, check_list)
    
    with open(header_filename, 'w') as f:
        f.write(header_code)
    with open(cpp_filename, 'w') as f:
        f.write(cpp_code)


