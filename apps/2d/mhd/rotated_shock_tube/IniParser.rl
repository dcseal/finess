#ifndef _INIREADER_H_
#define _INIREADER_H_
#include <iostream>
#include <streambuf>
#include <string>
#include <map>

#include "util.h"

%%{
    machine ini_parser;
    write data;

    action handle_error{
	fbreak;
    }


    action enter_section_name{
        p_important_token_start = p;
    }

    action exit_section_name{
        current_section_name = get_important_token();
    }

    action enter_item_name{
        p_important_token_start = p;
    }

    action exit_item_name{
        item_name = get_important_token();
    }

    action enter_item_content_plus_trailing_spaces{
        p_important_token_start = p;
    }

    action exit_item_content_plus_trailing_spaces{
        {
	    std::string item_content_plus_trailing_spaces =
	        get_important_token();
            item_content =
	        remove_trailing_spaces(item_content_plus_trailing_spaces);
	}
    }

    action exit_item_line{
        ini_doc[current_section_name][item_name] = item_content; 
    }

    action enter_line{
        ++line_number;
    } 

    non_newline_space = [ \t];
    spaces = non_newline_space*;
    newline = '\n' | '\r';
    comment_marker = ';' ;
    comments = comment_marker (any - newline)*;
    comments_line = spaces comments;

    section_name = (alpha (alnum | '_')*) 
                   >enter_section_name
		   %exit_section_name;
    section_line = spaces '[' section_name ']' spaces comments?;

    item_name = (alpha ( alnum | '-' | '_' | '(' | ')' )*)
                >enter_item_name
		%exit_item_name;
    item_content_plus_trailing_spaces = 
        (any - (comment_marker | newline) )*
	>enter_item_content_plus_trailing_spaces
	%exit_item_content_plus_trailing_spaces;
    item_line = (spaces item_name spaces '=' spaces 
                 item_content_plus_trailing_spaces
		 comments?)
		%exit_item_line;
    line =  (comments_line |
             section_line |
	     item_line |
	     '' )
	     >enter_line;


    main := (line newline)* (line newline?)?
            $!handle_error;
}%%

///@brief A simple .ini file parser class.
///
/// Example:
/// a.ini
/// @code{.ini}
///   [section_1]   ;some comments
///  a(234)_-33 = kk ll    ;some comments
///  abc = dd
///
///  [section_2]
///    bb = 2.3e9  ; some comments
///      [section_1]
///  a = dd
///  abc = ee
///  d = 
/// @endcode
///
/// C++ code
/// @code{.cpp}
/// IniParser parser;
/// std::string a_ini_file_content;
/// //read entire a.ini into a_ini_file_content
/// // ...
/// int parse_result = parser.parse(a_ini_file_content);
/// assert(parse_result == 0);
///
/// IniParser::ini_doc_type& ini_doc = parser.get_ini_doc_for_overwrite();
/// ini_doc["section_1"]["bcd"] = "efg";
///
/// assert(ini_doc["section_1"]["a(234)_-33"] == "kk ll");
/// assert(ini_doc["section_1"]["abc"] == "ee");
/// assert(ini_doc["section_1"]["a"] == "dd");
/// assert(ini_doc["section_1"]["d"] == "");
/// assert(ini_doc["section_1"]["bcd"] == "efg");
/// assert(ini_doc["section_2"]["bb"] == "2.3e9";
/// 
/// // parser.ini_doc_string() will return:
/// // [abc]
/// // bcd=eff
/// // [section_1]
/// // a=dd
/// // a(234)_-33=kk ll
/// // abc=ee
/// // d=
/// // [section_2]
/// // bb=2.3e9
/// @endcode

class IniParser{

private:
    const char *p, *pe;
    const char *eof;
    int cs;

    const char *p_important_token_start;

    int line_number;

    std::string current_section_name;
    std::string item_name;
    std::string item_content;
private:
    std::string get_important_token(){
        return std::string(p_important_token_start,
	                   p - p_important_token_start);
    }
private:

public:    
    typedef std::map<std::string, std::string> ini_section_type;
    typedef std::map<std::string, ini_section_type> ini_doc_type;

private:
    ini_doc_type ini_doc;

public:
    ///@brief Returns a read-only reference for the internal ini_doc
    ///       variable.
    ///
    ///@note Due to const constraint on the [] operator of std::map,
    ///      the caller cannot use a const ini_doc_type& variable to
    ///      hold the return value.  The caller must copy (possibly
    ///      implicitly via assignment) the return value to a local
    ///      variable in order to call [] operator later.
    ///      If speed is a concern, use get_ini_doc_for_overwrite()
    ///      instead.
    const ini_doc_type& get_ini_doc(){
        return ini_doc;
    }

    ///@brief Returns a writable reference for the internal ini_doc
    /// variable.
    ini_doc_type& get_ini_doc_for_overwrite(){
        return ini_doc;
    }

    ///@brief Returns a string-representation of the internal ini_doc
    /// variable.
    std::string ini_doc_as_string(){
        std::string ret;
        for(ini_doc_type::const_iterator it1 = ini_doc.begin();
            it1 != ini_doc.end();
           ++it1){
    	   ret.append("[");
	   ret.append(it1->first);
	   ret.append("]\n");
    	   for(ini_section_type::const_iterator it2 = (it1->second).begin();
    	       it2 != (it1->second).end();
    	       ++it2){
    	       ret.append(it2->first);
	       ret.append("=");
	       ret.append(it2->second);
	       ret.append("\n");
    	   }
        } 
	return ret;
    }

public:   
    ///@brief Parses the content of a .ini file.
    ///
    /// get_ini_doc_for_overwrite(), or ini_doc_as_string()
    /// To get the result of parsing, use get_ini_doc(),
    ///@return 0 if parsing is successful; 
    ///        otherwise the line number of first error.
    int parse(const std::string& ini_file_content){
        using std::string;
        p = ini_file_content.c_str();
        pe = p + ini_file_content.length();
        eof = 0;
    
        line_number = 0;
        // types of important_token:
        // - section name
        // - item name
        // - item content (with trailing spaces)
    
        %%write init;
        %%write exec;
        
        return cs >= ini_parser_first_final ? 0 : line_number;
    }

};


#endif
