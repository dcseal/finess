#ifndef _TESTINIPARSER_H_
#define _TESTINIPARSER_H_

#include <stdexcept>
#include <iostream>
#include <fstream>

#include <cxxtest/TestSuite.h>

#include "IniParser.h"

class TestIniParser: public CxxTest::TestSuite{
public:
    class Setup{
	public:
	    std::string filename;
	    Setup(const std::string& filename, const std::string& content)
	    {	
		this->filename = filename;
	        while(existFile(this->filename)){
	            this->filename += "svuf";
	        }
		CxxTest::setAbortTestOnFail(true);
		TSM_ASSERT(filename 
			   + " expected to be non-existent! "
			     " Unable to continue tests. "
			     " Possible solutions: remove the file, "
			     "or change its name",
			   !existFile(filename));
		CxxTest::setAbortTestOnFail(CXXTEST_DEFAULT_ABORT);
		          
		std::ofstream ofs(this->filename.c_str());
		ofs << content;
		ofs.close();
	    }
	    ~Setup(){
		std::remove(this->filename.c_str());		    
	    }
    };

    void testIniParser(){
	{
            std::string content;
	    content += "[section_1]   ;some comments\n";
	    content += "a(234)_-33 = kk ll    ;some comments\n";
	    content += "abc = dd\n";
	    content += "[section_2]\n";
	    content += "\n";
	    content += "  bb = 2.3e9  ; some comments\n";
	    content += "    [section_1]\n";
	    content += "a = dd\n";
	    content += "abc = ee\n";
	    content += "d = ;\n";
            
	    Setup setup("someveryunlikelyfilename11223344", content);
	    
	    std::string filename = setup.filename;
            std::string ini_file_content;
	    {
		std::ifstream ifs(filename.c_str());
		ini_file_content = read_entire_stream(ifs);
	    }
	    IniParser parser;
	    parser.parse(ini_file_content);
	    IniParser::ini_doc_type& ini_doc_reference = 
		parser.get_ini_doc_for_overwrite();
	    ini_doc_reference["section_1"]["bcd"] = "efg";

            IniParser::ini_doc_type ini_doc = parser.get_ini_doc();

            TSM_ASSERT(content, ini_doc["section_1"]["a(234)_-33"] == "kk ll");
            TSM_ASSERT(content, ini_doc["section_1"]["abc"] == "ee");
            TSM_ASSERT(content, ini_doc["section_1"]["a"] == "dd");
            TSM_ASSERT(content, ini_doc["section_1"]["d"] == "");
            TSM_ASSERT(content, ini_doc["section_1"]["bcd"] == "efg");
            TSM_ASSERT(content, ini_doc["section_2"]["bb"] == "2.3e9");
	}
	{
	    IniParser::ini_doc_type ini_doc_input;
	    ini_doc_input["section_1"]["params1"] = "123";
	    ini_doc_input["section_2"]["params2"] = "u-7";
	    std::string content = IniParser::ini_doc_as_string(ini_doc_input);

	    Setup setup("someveryunlikelyfilename11223344", content);
            
	    std::string filename = setup.filename;
            std::string ini_file_content;
	    {
		std::ifstream ifs(filename.c_str());
		ini_file_content = read_entire_stream(ifs);
	    }
	    IniParser parser;
	    parser.parse(ini_file_content);
            IniParser::ini_doc_type ini_doc = parser.get_ini_doc();

	    TSM_ASSERT(content, ini_doc == ini_doc_input);


	}
    }

};

#endif
