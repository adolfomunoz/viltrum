#include <iostream>
#include <iomanip>
#include "../quadrature/integrate.h"
#include "../quadrature/integrate-bins.h"
#include "../quadrature/integrate-bins-adaptive.h"
#include "../quadrature/integrate-optimized-adaptive-stratified-control-variates.h"
#include "../quadrature/monte-carlo.h"
#include <cmath>
#include "../functions/functions1d.h"
//#include <numeric>
//#include <execution>

template<typename F>
class FunctionWrapper {
	F f;
	mutable unsigned long evals = 0;
public:
	FunctionWrapper(const F& f) : f(f) {}
	FunctionWrapper(F&& f) : f(std::forward<F>(f)) {}
	
    template<typename Float>
	auto operator()(const std::array<Float, 1>& pos) const {
		++evals;
		return f(std::get<0>(pos));
	}
	
	unsigned long evaluations() const { return evals; }
};

template<typename F, typename I>
void test(const char* name, const I& integrator, unsigned int bins, const F& f, double rmin = 0, double rmax=1) {
	auto wrapper = FunctionWrapper<F>(f);
    std::vector<double> result(bins);
    integrate_bins(integrator,result,wrapper, range(rmin, rmax));
	std::cout<<name<<" -> [";
    for (double r : result)
        std::cout<<std::setw(10)<<std::setprecision(2)<<std::scientific<<r<<" ";
    std::cout <<"]\t("<<std::setw(6)<<wrapper.evaluations()<<" evaluations)"<<std::endl;
}

int main(int argc, char **argv) {
	auto [f,groundtruth] = function1d(argc,argv);
    unsigned int bins = 4;
	double rmin = 0; double rmax = 1;
	for (int i = 0; i < argc; ++i) {
		if ( (std::string(argv[i]) == "-range") && (i<(argc-2)) ) {
			rmin = atof(argv[++i]); rmax = atof(argv[++i]);
		}
        else if ( (std::string(argv[i]) == "-bins") && (i<(argc-1)) ) {
            bins = atoi(argv[++i]);
        }
	}

	std::cout<<"Ground truth  "<<" -> [";
    double dr = (rmax - rmin)/double(bins);
    for (unsigned int i = 0; i < bins; ++i)
        std::cout<<std::setw(10)<<std::setprecision(2)<<std::scientific<<bins*groundtruth(rmin + i*dr,rmin + (i+1)*dr)<<" ";
    std::cout<<" ]"<<std::endl;
	std::cout<<std::endl<<"Per pixel"<<std::endl;	
	test("Trapezoidal   ",integrator_bins_per_bin(integrator_quadrature(trapezoidal)), bins, f,rmin,rmax);
	test("Simpson       ",integrator_bins_per_bin(integrator_quadrature(simpson)), bins, f,rmin,rmax);
	test("Boole         ",integrator_bins_per_bin(integrator_quadrature(boole)), bins, f,rmin,rmax);
	test("M-C         1 ",integrator_bins_per_bin(integrator_monte_carlo_uniform(1)), bins, f,rmin,rmax);
	test("M-C        10 ",integrator_bins_per_bin(integrator_monte_carlo_uniform(10)), bins, f,rmin,rmax);
	test("M-C       100 ",integrator_bins_per_bin(integrator_monte_carlo_uniform(100)), bins, f,rmin,rmax);
	test("M-C      1000 ",integrator_bins_per_bin(integrator_monte_carlo_uniform(1000)), bins, f,rmin,rmax);
	std::cout<<std::endl<<"Global"<<std::endl;	
	test("Simpson       ",integrator_bins_stepper(stepper_bins_adaptive(nested(simpson,trapezoidal)),0), bins, f,rmin,rmax);
	test("Boole         ",integrator_bins_stepper(stepper_bins_adaptive(nested(boole,simpson)),0), bins, f,rmin,rmax);
	test("M-C         1 ",integrator_bins_stepper(stepper_bins_monte_carlo_uniform(),1), bins, f,rmin,rmax);
	test("M-C        10 ",integrator_bins_stepper(stepper_bins_monte_carlo_uniform(),10), bins, f,rmin,rmax);
	test("M-C       100 ",integrator_bins_stepper(stepper_bins_monte_carlo_uniform(),100), bins, f,rmin,rmax);
	test("M-C      1000 ",integrator_bins_stepper(stepper_bins_monte_carlo_uniform(),1000), bins, f,rmin,rmax);
	
	std::cout<<std::endl<<"Adaptive optimized stratified control variates"<<std::endl;
	test("ST   1      4 ",integrator_optimized_adaptive_stratified_control_variates(
		nested(simpson,trapezoidal),error_relative_single_dimension_size(1.e-6),1,4), bins, f, rmin, rmax);
	test("ST   2      4 ",integrator_optimized_adaptive_stratified_control_variates(
		nested(simpson,trapezoidal),error_relative_single_dimension_size(1.e-6),2,4), bins, f, rmin, rmax);
	test("ST   4      4 ",integrator_optimized_adaptive_stratified_control_variates(
		nested(simpson,trapezoidal),error_relative_single_dimension_size(1.e-6),4,4), bins, f, rmin, rmax);
	test("ST   4      8 ",integrator_optimized_adaptive_stratified_control_variates(
		nested(simpson,trapezoidal),error_relative_single_dimension_size(1.e-6),4,8), bins, f, rmin, rmax);
	test("ST   8      8 ",integrator_optimized_adaptive_stratified_control_variates(
		nested(simpson,trapezoidal),error_relative_single_dimension_size(1.e-6),8,8), bins, f, rmin, rmax);

	return 0;
}



