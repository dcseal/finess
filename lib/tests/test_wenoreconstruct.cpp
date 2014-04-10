#include <assert.h>
#include <vector>
#include <utility>
#include <iostream>
#include <cmath>
#include <string>

#include <functional>

#include "tensors.h"
#include "WenoReconstruct.h"
#include "constants.h"

//typedef std::function<void(const dTensor2&, dTensor2&)> reconstruct_t;
//typedef void (*reconstruct_t)(const dTensor2&, dTensor2&);

//typedef std::function<double(double)> realfunc_t;
typedef double (*realfunc_t)(double);

///@brief Computes derivative of f at x, using reconstruction R, assuming using the 'plus' version.
double derivative(realfunc_t f, double x,
				  reconstruct_t R, int order, double dx){
	using std::vector;
	assert(order%2 == 1);
	int r = (order + 1)/2;
	vector<double> sample_f(2*r);
	for(int i = 0; i < 2*r; ++i){
		sample_f[i] = f(x + (i - r) * dx);
	}
	
	dTensor2 g(1, order);
	dTensor2 diff_g(1, 1);
	double fm, fp;
	
	for(int i = 0; i < order; ++i){
		g.set(1, i+1, sample_f[i]);
	}
	R(g, diff_g);
	fm = diff_g.get(1, 1);
	
	for(int i = 0; i < order; ++i){
		g.set(1, i+1, sample_f[i+1]);
	}
	R(g, diff_g);
	fp = diff_g.get(1, 1);
	//std::cout << fm << std::endl;
	//std::cout << fp << std::endl;
	return (fp - fm) / dx;
}

double derivative_error(realfunc_t f, realfunc_t df, double x, reconstruct_t R, int order, double dx){
	double df_approx = derivative(f, x, R, order, dx);
	double df_exact  = df(x);
	return std::abs(df_approx - df_exact);
}


struct WenoReconstructTestCase{
	std::string  docString; ///<@brief Description of the test case.  Will be included in output.
	reconstruct_t R;  ///<@brief Reconstruction function to be tested
	int order;        ///<@brief Order of the reconstruction method.
					  ///  Used to determine the number of stencil points.
	realfunc_t f;     ///<@brief A real to real function.
	realfunc_t df;	  ///<@brief The exact derivative function of <tt>f</tt>.
	double x;		  ///<@brief The real number at which the estimated and exact derivatives are compared.
	double init_dx;   ///<@brief Initial @f$@\Delta xf$.
	double shrinking_factor;  ///<@brief The factor by which @f$@\Delta xf$ gets divided each time.
	int nRuns;		  ///<@brief Number of runs.
};


double f(double x){
	return sin(50*x);
}

double df(double x){
	return cos(50*x)*50;
}


double g(double x){
	return exp(pow(1-sin(1/x), 2));
}

double dg(double x){
	return 2*exp(pow(1-sin(1/x), 2)) * cos(1/x) * (1-sin(1/x)) / (x*x);
}


///@brief The test cases.
WenoReconstructTestCase testCases[]
={
	{
		"JS 5, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_JS5,
		5,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"JS 7, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_JS7,
		7,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"JS 9, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_JS9,
		9,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"JS 11, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_JS11,
		11,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"Z 5, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_Z5,
		5,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"Z 7, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_Z7,
		7,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"Z 9, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_Z9,
		9,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"FD 5, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_FD5,
		5,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"FD 7, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_FD7,
		7,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},
	{
		"FD 9, sin(50x), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_FD9,
		9,
		f,
		df,
		pi/100,
		1.0/8.0,
		2.0,
		20
	},

///////////////////////////////////////////////////////////////////////////////
//  A different function, and different x.

	{
		"JS 5, exp((1-sin(1/x))^2)), at pi/100, init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_JS5,
		5,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"JS 7, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_JS7,
		7,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"JS 9, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_JS9,
		9,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"JS 11, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_JS11,
		11,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"Z 5, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_Z5,
		5,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"Z 7, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_Z7,
		7,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"Z 9, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_Z9,
		9,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"FD 5, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_FD5,
		5,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"FD 7, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_FD7,
		7,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	},
	{
		"FD 9, exp((1-sin(1/x))^2)), at 1/(10.5*pi), init_dx=1/8, shrinking_factor=2, nRuns=20",
		WenoReconstruct_FD9,
		9,
		g,
		dg,
		1/(10.5*pi),
		1.0/8.0,
		2.0,
		20
	}
};


void test(WenoReconstructTestCase testCase){
	using std::log;

	// double x = pi/100;
	// reconstruct_t R = WenoReconstruct_JS9;
	// int order = 9;

	// double init_dx = 1.0 / 8.0;
	// double shrinking_factor = 2.0;
	// assert(shrinking_factor > 1.0);
	// int nRuns = 20;

	// double dx = init_dx;

	std::cout << "BEGIN --- " << testCase.docString << std::endl;
	std::cout << "-------------------------" << std::endl;

	double x = testCase.x;
	reconstruct_t R = testCase.R;
	int order = testCase.order;
	double init_dx = testCase.init_dx;
	double shrinking_factor = testCase.shrinking_factor;
	assert(shrinking_factor > 1.0);
	int nRuns = testCase.nRuns;
	realfunc_t f = testCase.f;
	realfunc_t df = testCase.df;

	double dx = init_dx;
	double error_old = NAN, error_new = NAN;
	for(int i = 0; i < nRuns; ++i){
		error_old = error_new;
		error_new = derivative_error(f, df, x, R, order, dx);
		std::cout << "dx = " << dx << std::endl;
		std::cout << "Old error = " << error_old << std::endl;
		std::cout << "New error = " << error_new << std::endl;
		std::cout << "Order = " << log(error_old / error_new) / log(shrinking_factor) << std::endl;
		std::cout << "-------------------------" << std::endl;
		dx /= shrinking_factor;
	}
	std::cout << "END   --- " << testCase.docString << std::endl;
	std::cout << "**********************************************************************************" << std::endl;
}


int main(){
	const int nTestCases = sizeof(testCases) / sizeof(WenoReconstructTestCase);
	for(int i = 0; i < nTestCases; ++i){
		test(testCases[i]);
	}
}