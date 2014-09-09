#ifndef _TESTPARAMS_H_
#define _TESTPARAMS_H_

#include <cxxtest/TestSuite.h>
#include "Params.h"

class TestParams: public CxxTest::TestSuite{
public:
    class Setup{
	public:
	    const string filename;
	    IniDocument ini_doc;
	    Setup(const string& filename, const string& content):
		filename(filename)
	    {
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
	while(existFile(filename)){
	    filename += "svuf";
	}
	CxxTest::setAbortTestOnFail(true);
	TSM_ASSERT(filename 
		   + " expected to be non-existent! "
		     " Unable to continue tests. "
		     " Possible solutions: remove the file, "
		     "or change its name",
		   !existFile(filename));
	CxxTest::setAbortTestOnFail(CXXTEST_DEFAULT_ABORT);


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
	
	//Parsing error will abort the program.  Omit test.
    }

};

#endif
