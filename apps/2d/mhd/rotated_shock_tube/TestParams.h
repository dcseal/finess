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
	    }
	    ~Setup(){
		std::remove(this->filename.c_str());		    
	    }
    };
	
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
	    content += " = \n";
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
		    params.get_meqn(),
		    2);

	    TSM_ASSERT_EQUALS(content, 
		    params.get_meqn_times_2(),
		    4);

	    TSM_ASSERT_EQUALS(content,
		    params.get_default_1(),
		    1);
	
	}

	{
	    string content;
	    content += "[reconstruction]\n";
	    content += "method = B\n";
	    content += "[dogParams]\n";
	    content += "meqn = 2 \n";
	    content += "one_of_1_3_5 = 1\n";
	    content += "default_1 = 2\n";
	    
	    Setup setup(filename, content);
	    Params params;

	    params.init(setup.filename);

	    TSM_ASSERT_EQUALS(content,
		    params.get_meqn(),
		    2);

	    TSM_ASSERT_EQUALS(content, 
		    params.get_meqn_times_2(),
		    4);

	    TSM_ASSERT_EQUALS(content,
		    params.get_default_1(),
		    2);
	
	}
    }

};

#endif
