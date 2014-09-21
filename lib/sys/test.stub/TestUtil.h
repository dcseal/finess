#ifndef _TESTUTIL_H_
#define _TESTUTIL_H_

#include <cxxtest/TestSuite.h>

#include <fstream>

#include "util.h"

class TestUtil: public CxxTest::TestSuite{
public:
    void testStringToAny(){
	using std::string;
	TS_ASSERT(stringToAny<string>("astring") == "astring");
	TS_ASSERT(stringToAny<string>("  astring") == "astring");
	TS_ASSERT(stringToAny<string>("  astring  ") == "astring");
	TS_ASSERT(stringToAny<string>("\"\"") == "\"\"");
	TS_ASSERT(stringToAny<string>("\"a b\"") == "\"a");

	TS_ASSERT(stringToAny<int>("12") == 12);
	TS_ASSERT(stringToAny<int>("+12") == 12);
	TS_ASSERT(stringToAny<int>("-12") == -12);

	TS_ASSERT(stringToAny<int>("1.6") == 1);
	TS_ASSERT(stringToAny<int>("-2.5") == -2);
	TS_ASSERT(stringToAny<int>("abc") == 0);
	TS_ASSERT(stringToAny<int>("++--**//") == 0);

	TS_ASSERT(stringToAny<int>("  12") == 12);
	TS_ASSERT(stringToAny<int>("  12abc") == 12);

	const double tolerance = 1e-12;
	TS_ASSERT_DELTA(stringToAny<double>("0.0"), 0.0, tolerance);
	TS_ASSERT_DELTA(stringToAny<double>("-0.0"), -0.0, tolerance);
	TS_ASSERT_DELTA(stringToAny<double>("1.2"), 1.2, tolerance);
	TS_ASSERT_DELTA(stringToAny<double>("+1.2"), 1.2, tolerance);
	TS_ASSERT_DELTA(stringToAny<double>("-1.2"), -1.2, tolerance);

	TS_ASSERT_DELTA(stringToAny<double>("1e-8"), 1e-8, tolerance);
	TS_ASSERT_DELTA(stringToAny<double>("-1.2e2"), -120, tolerance);

    }

    void testAnyToString(){
	using std::string;

	TS_ASSERT(anyToString(123) == "123");
	TS_ASSERT(anyToString(-123) == "-123");
    }

    void testExistFile(){
	using std::string;
	using std::ofstream;
	const string unlikely_filename = "3844someveryunlikelyfilename8922";
    	class Setup{
	    public:
		const string filename;
		Setup(const string& filename):
		    filename(filename)
	        {
		    using std::ofstream;
		    ofstream(this->filename.c_str());
		}
		~Setup(){
		    using std::remove;
		    remove(this->filename.c_str());
		}
	};
	
	string filename = unlikely_filename;

	CxxTest::setAbortTestOnFail(true);
	TSM_ASSERT(filename 
		   + " expected to be non-existent! "
		     " Unable to continue tests. "
		     " Possible solutions: remove the file, "
		     "or change its name",
		   !existFile(filename));
	CxxTest::setAbortTestOnFail(CXXTEST_DEFAULT_ABORT);

	{Setup setup(filename);	    	    	
	    TS_ASSERT(existFile(filename));
	}

	TS_ASSERT(!existFile(filename));
    }

};


#endif
