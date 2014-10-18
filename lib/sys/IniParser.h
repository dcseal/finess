
#line 1 "IniParser.rl"
#ifndef _INIREADER_H_
#define _INIREADER_H_
#include <iostream>
#include <streambuf>
#include <string>
#include <map>

#include "util.h"


#line 14 "IniParser.h"
static const char _ini_parser_actions[] = {
	0, 1, 0, 1, 1, 1, 2, 1, 
	3, 1, 4, 1, 5, 1, 6, 1, 
	7, 1, 8, 2, 5, 6, 2, 6, 
	7, 2, 8, 3, 3, 5, 6, 7
	
};

static const char _ini_parser_key_offsets[] = {
	0, 0, 8, 21, 24, 28, 36, 46, 
	48, 53, 56, 58
};

static const char _ini_parser_trans_keys[] = {
	9, 32, 59, 91, 65, 90, 97, 122, 
	9, 32, 45, 61, 95, 40, 41, 48, 
	57, 65, 90, 97, 122, 9, 32, 61, 
	65, 90, 97, 122, 93, 95, 48, 57, 
	65, 90, 97, 122, 9, 10, 13, 32, 
	59, 91, 65, 90, 97, 122, 10, 13, 
	9, 10, 13, 32, 59, 10, 13, 59, 
	10, 13, 9, 10, 13, 32, 59, 0
};

static const char _ini_parser_single_lengths[] = {
	0, 4, 5, 3, 0, 2, 6, 2, 
	5, 3, 2, 5
};

static const char _ini_parser_range_lengths[] = {
	0, 2, 4, 0, 2, 3, 2, 0, 
	0, 0, 0, 0
};

static const char _ini_parser_index_offsets[] = {
	0, 0, 7, 17, 21, 24, 30, 39, 
	42, 48, 52, 55
};

static const char _ini_parser_indicies[] = {
	1, 1, 2, 4, 3, 3, 0, 5, 
	5, 6, 7, 6, 6, 6, 6, 6, 
	0, 8, 8, 9, 0, 10, 10, 0, 
	12, 11, 11, 11, 11, 0, 13, 14, 
	14, 13, 15, 17, 16, 16, 0, 18, 
	18, 2, 20, 21, 21, 20, 22, 19, 
	24, 24, 25, 23, 27, 27, 26, 28, 
	18, 18, 28, 2, 0, 0
};

static const char _ini_parser_trans_targs[] = {
	0, 1, 7, 2, 4, 3, 2, 8, 
	3, 8, 5, 5, 11, 1, 6, 7, 
	2, 4, 6, 9, 8, 6, 10, 9, 
	6, 10, 10, 6, 11
};

static const char _ini_parser_trans_actions[] = {
	1, 0, 0, 7, 0, 9, 0, 9, 
	0, 0, 3, 0, 5, 17, 17, 17, 
	25, 17, 0, 11, 11, 28, 19, 0, 
	22, 13, 0, 15, 0
};

static const char _ini_parser_eof_actions[] = {
	0, 1, 1, 1, 1, 1, 17, 0, 
	28, 22, 15, 0
};

static const int ini_parser_start = 6;
static const int ini_parser_first_final = 6;
static const int ini_parser_error = 0;

static const int ini_parser_en_main = 6;


#line 88 "IniParser.rl"


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

    ///@brief Returns a string-representation of an ini_doc variable. 
    static std::string ini_doc_as_string(ini_doc_type& ini_doc){
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
    
        
#line 232 "IniParser.h"
	{
	cs = ini_parser_start;
	}

#line 228 "IniParser.rl"
        
#line 239 "IniParser.h"
	{
	int _klen;
	unsigned int _trans;
	const char *_acts;
	unsigned int _nacts;
	const char *_keys;

	if ( p == pe )
		goto _test_eof;
	if ( cs == 0 )
		goto _out;
_resume:
	_keys = _ini_parser_trans_keys + _ini_parser_key_offsets[cs];
	_trans = _ini_parser_index_offsets[cs];

	_klen = _ini_parser_single_lengths[cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += (unsigned int)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _ini_parser_range_lengths[cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += (unsigned int)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	_trans = _ini_parser_indicies[_trans];
	cs = _ini_parser_trans_targs[_trans];

	if ( _ini_parser_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _ini_parser_actions + _ini_parser_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 14 "IniParser.rl"
	{
	{p++; goto _out; }
    }
	break;
	case 1:
#line 19 "IniParser.rl"
	{
        p_important_token_start = p;
    }
	break;
	case 2:
#line 23 "IniParser.rl"
	{
        current_section_name = get_important_token();
    }
	break;
	case 3:
#line 27 "IniParser.rl"
	{
        p_important_token_start = p;
    }
	break;
	case 4:
#line 31 "IniParser.rl"
	{
        item_name = get_important_token();
    }
	break;
	case 5:
#line 35 "IniParser.rl"
	{
        p_important_token_start = p;
    }
	break;
	case 6:
#line 39 "IniParser.rl"
	{
        {
	    std::string item_content_plus_trailing_spaces =
	        get_important_token();
            item_content =
	        remove_trailing_spaces(item_content_plus_trailing_spaces);
	}
    }
	break;
	case 7:
#line 48 "IniParser.rl"
	{
        ini_doc[current_section_name][item_name] = item_content; 
    }
	break;
	case 8:
#line 52 "IniParser.rl"
	{
        ++line_number;
    }
	break;
#line 372 "IniParser.h"
		}
	}

_again:
	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	const char *__acts = _ini_parser_actions + _ini_parser_eof_actions[cs];
	unsigned int __nacts = (unsigned int) *__acts++;
	while ( __nacts-- > 0 ) {
		switch ( *__acts++ ) {
	case 0:
#line 14 "IniParser.rl"
	{
	{p++; goto _out; }
    }
	break;
	case 5:
#line 35 "IniParser.rl"
	{
        p_important_token_start = p;
    }
	break;
	case 6:
#line 39 "IniParser.rl"
	{
        {
	    std::string item_content_plus_trailing_spaces =
	        get_important_token();
            item_content =
	        remove_trailing_spaces(item_content_plus_trailing_spaces);
	}
    }
	break;
	case 7:
#line 48 "IniParser.rl"
	{
        ini_doc[current_section_name][item_name] = item_content; 
    }
	break;
	case 8:
#line 52 "IniParser.rl"
	{
        ++line_number;
    }
	break;
#line 423 "IniParser.h"
		}
	}
	}

	_out: {}
	}

#line 229 "IniParser.rl"
        
        return cs >= ini_parser_first_final ? 0 : line_number;
    }

};


#endif
