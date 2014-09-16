


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
         "default_value": str(default_value)}
    return default_on_missing



class ParameterType:
    def __init__(self, type_string):
        if type(type_string) != str:
            raise ValueError("type_string should be a string.")                    
                             
        canonical_types = {"int", "double", "string"}
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
        self.full_enum_scope_name = class_name + "::" + enum_scope_name
        self.string_enumerator_dict = string_enumerator_dict
        self.type_string = """%(class_name)s::%(enum_scope_name)s::enum_type""" %                            {"class_name": class_name,
                            "enum_scope_name": enum_scope_name}
        
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
    def __init__(self):
        pass
    
    def get_cpp_code(self):
        return ""


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


def generate_header_cpp(class_name, global_variable_name, parameters, accessors, checks):
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
    types_requiring_definitions =         [t for t in set([p.type_ for p in parameters])
             if getattr(t, "generate_type_definition", None) != None]
   
    class_declaration_contents =         "// Type definitions \n" +         "\n".join([t.generate_type_definition() for t in types_requiring_definitions]) +         "\n\n// Member variables declarations\n" +         "\n".join([p.get_declaration_code() for p in parameters]) +         "\n\n// Accessor definitions\n" +         "\n".join([a.get_accessor_code() for a in accessors])
        
        
    types_requiring_stringToAny_specialization =         [t for t in set([p.type_ for p in parameters])
             if getattr(t, "generate_stringToAny_specialization", None) != None]
    stringToAny_specializations_code =         "//stringToAny specializations\n" +         "\n".join([t.generate_stringToAny_specialization() for t in types_requiring_stringToAny_specialization])
    
    class_declaration = """
class Params{
public:
    Params(){}
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
{"class_declaration_contents": class_declaration_contents}

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
void Params::init(const std::string& inputFilename){
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
{"cpp_init_method_contents": cpp_init_method_contents}

    header_code = header_head + class_declaration + stringToAny_specializations_code + header_tail
    cpp_code = cpp_head + cpp_init_method
    
    return header_filename, header_code, cpp_filename, cpp_code
    

    

if __name__ == "__main__":
    class_name = "Params"
    global_variable_name = "params"
    
    
    parameters = [Parameter(variable_name = "reconstruction_method",
                            section = "reconstruction",
                            name = "method",
                            type_ = EnumParameterType(enum_scope_name = "ReconstructionMethod",
                                                      string_enumerator_dict = \
                                                          {"A": "A", "B":"B", "C":"C"})),
                  Parameter(variable_name = "meqn",
                            section = "dogParams",
                            name = "meqn",
                            type_ = "int"),
                  Parameter(variable_name = "one_of_1_3_5", 
                            section = "dogParams",
                            name = "one_of_1_3_5", 
                            type_ = "int"),
                  DerivedParameter(variable_name = "meqn_times_2",
                                   type_ = "int",
                                   defining_expression_in_cpp = "this->meqn * 2"),
		  Parameter(variable_name = "default_1",
		            section = "dogParams",
			    name = "default_1",
			    type_ = "int",
			    default_value = 1)]
    
    accessors = map(Accessor, parameters)
    
    checks = [CheckGreaterThan(parameters[1], 0),
              CheckOneOf(parameters[2], [1, 3, 5])]
    
    header_filename, header_code, cpp_filename, cpp_code = \
        generate_header_cpp(class_name, global_variable_name, parameters, accessors, checks)
    
    with open(header_filename, 'w') as f:
        f.write(header_code)
    with open(cpp_filename, 'w') as f:
        f.write(cpp_code)
    
    print header_code
    
    print cpp_code
