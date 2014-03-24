#include <assert.h>
#include <vector>
#include <utility>
#include <iostream>
#include <cmath>

#include <functional>

#include "tensors.h"
#include "WenoReconstruct.h"


//typedef std::function<void(const dTensor2&, dTensor2&)> reconstruct_t;
typedef void (*reconstruct_t)(const dTensor2&, dTensor2&);

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

double f(double x){
	return sin(50*x);
}

double df(double x){
	return cos(50*x)*50;
}

void test(){
	using std::log;

	//realfunc_t f = [](double x){return sin(50*x);}, df = [](double x){return cos(50*x)*50;};
	double x = 5.0;
	reconstruct_t R = WenoReconstruct_JS7;
	int order = 7;

	double init_dx = 1.0 / 8.0;
	double shrinking_factor = 2.0;
	assert(shrinking_factor > 1.0);
	int nRuns = 20;

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
}


int main(){
	test();	
}