#ifndef _TESTPARAMS_H_
#define _TESTPARAMS_H_

#include <stdexcept>

#include <cxxtest/TestSuite.h>
#include "Params.h"

class TestParams: public CxxTest::TestSuite{
public:
    class Setup{
	public:
	    string filename;
	    IniDocument ini_doc;
	    Setup(const string& filename, const string& content)
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
		          
		ofstream ofs(this->filename.c_str());
		ofs << content;
		ofs.close();
		ini_doc.initFromFile(this->filename);
	    }
	    ~Setup(){
		std::remove(this->filename.c_str());		    
	    }
    };
	
    void testIniDocument(){
	using std::string;
	using std::ofstream;

	string filename = "3844someveryunlikelyfilename8922";
	{
	    string content;
	    content += "[section1]\n";
	    content += "a = b\n";
	    content += "c = d\n";
	    Setup setup(filename, content);
	    IniDocument& ini_doc = setup.ini_doc;
	    
	    TSM_ASSERT(content, ini_doc["section1"]["a"] == "b");
	    TSM_ASSERT(content, ini_doc["section1"]["c"] == "d");
	}

	{
	    string content;
	    content += "[section1]\n";
	    content += "a(b)  = 1e08    ;some comments\n";
	    content += "a(b)  = 1e09\n";
	    content += "\n";
	    content += ";entire line of comments\n";
	    content += "[section2]\n";
	    content += "a(b)  = 12  ;something else\n";

	    Setup setup(filename, content);
	    IniDocument& ini_doc = setup.ini_doc;
	
	    TSM_ASSERT_EQUALS(content,
		    ini_doc["section1"]["a(b)"],
		    "1e08");
	    //Note the repeated entry was ignored.
	    TSM_ASSERT_EQUALS(content,
	            ini_doc["section2"]["a(b)"],
		    "12");
	    TSM_ASSERT_EQUALS(content, 
		    ini_doc["section3"]["a(b)"],
		    "");
	    TSM_ASSERT_EQUALS(content,
		    ini_doc["section2"]["nonexist"],
		    "");
	}

	{
	    string content;
	    content += "[section1]\n";
	    content += "a = b\n";
	    content += "c = d\n";
	    content += "[section2]\n";
	    content += "e = f\n";
	    content += "[section1]\n";
	    content += "a = g\n";
	    content += "h = i\n";

	    Setup setup(filename, content);
	    IniDocument& ini_doc = setup.ini_doc;

	    TSM_ASSERT_EQUALS(content,
		    ini_doc["section1"]["a"],
		    "b");
	    TSM_ASSERT_EQUALS(content,
		    ini_doc["section1"]["c"],
		    "d");
	    TSM_ASSERT_EQUALS(content,
		    ini_doc["section1"]["h"],
		    "");
	    TSM_ASSERT_EQUALS(content,
		    ini_doc["section2"]["e"],
		    "f");
	}
	
	//Parsing error will abort the program.  Omit test.
    }

    void testParams(){
	using std::string;


	string filename = "3844someveryunlikelyfilename8922";

	{
	    string content;
	    content += "[reconstruction]\n";
	    content += "method = B\n";
	    content += "[dogParams]\n";
	    content += "meqn = 2 \n";
	    content += "one_of_1_3_5 = 1\n";

	    Setup setup(filename, content);
	    Params params;
	    params.init(setup.filename);
	    
	    TSM_ASSERT_EQUALS(content, 
		    params.get_reconstruction_method(),
		    Params::ReconstructionMethod::B);
	    TSM_ASSERT_EQUALS(content,
		    params.get_meqn(),
		    2);	    
	}

	{
	    string content;
	    content += "[reconstruction]\n";
	    content += "method = someinvalidvalue\n";
	    content += "[dogParams]\n";
	    content += "meqn = 2 \n";
	    content += "one_of_1_3_5 = 1\n";

	    Setup setup(filename, content);
	    Params params;
	    TSM_ASSERT_THROWS(content,
		    params.init(setup.filename),
		    std::runtime_error);
	    
	}
	{
	    string content;
	    content += "[reconstruction]\n";
	    content += "method = A\n";
	    content += "[dogParams]\n";
	    content += "one_of_1_3_5 = 1\n";
            content += "\n";

	    Setup setup(filename, content);
	    Params params;
	    TSM_ASSERT_THROWS(content,
		    params.init(setup.filename),
		    std::runtime_error);
	}
	{
	    string content;
	    content += "[reconstruction]\n";
	    content += "method = A\n";
	    content += "[dogParams]\n";
	    content += "meqn = -1\n";
	    content += "one_of_1_3_5 = 1\n";
		
	    Setup setup(filename, content);
	    Params params;
	    TSM_ASSERT_THROWS(content,
		    params.init(setup.filename),
		    std::runtime_error);
    
	}
	{
	    string content;
	    content += "[reconstruction]\n";
	    content += "method = B\n";
	    content += "[dogParams]\n";
	    content += "meqn = 2 \n";
	    content += "one_of_1_3_5 = 1\n";
	    
	    Setup setup(filename, content);
	    Params params;

	    params.init(setup.filename);
	    TSM_ASSERT_EQUALS(content, 
		    params.get_one_of_1_3_5(),
		    1);
	}
	{
	    string content;
	    content += "[reconstruction]\n";
	    content += "method = B\n";
	    content += "[dogParams]\n";
	    content += "meqn = 2 \n";
	    content += "one_of_1_3_5 = 2\n";
	    
	    Setup setup(filename, content);
	    Params params;
	    TSM_ASSERT_THROWS(content, 
         	    params.init(setup.filename),
		    std::runtime_error);
	}
    }

};

#endif
