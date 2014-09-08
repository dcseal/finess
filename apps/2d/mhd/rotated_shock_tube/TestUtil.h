#ifndef _TESTUTIL_H_
#define _TESTUTIL_H_

#include <cxxtest/TestSuite.h>
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

};


#endif
