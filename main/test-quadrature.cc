#include <iostream>
#include <iomanip>
#include "../quadrature/integrate.h"
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
void test(const char* name, const I& integrator, const F& f, double rmin = 0, double rmax=1) {
	auto wrapper = FunctionWrapper<F>(f);
    auto integral = integrator.integrate(wrapper, range(rmin, rmax));
	std::cout<<name<<" -> "
	         <<std::setw(15)<<std::setprecision(8)<<std::scientific<<integral
			 <<"\t("<<std::setw(6)<<wrapper.evaluations()<<" evaluations)"<<std::endl;
}

int main(int argc, char **argv) {
	auto [f,groundtruth] = function1d(argc,argv);
	double rmin = 0; double rmax = 1;
	for (int i = 0; i < argc; ++i) {
		if ( (std::string(argv[i]) == "-range") && (i<(argc-2)) ) {
			rmin = atof(argv[++i]); rmax = atof(argv[++i]);
		}
	}


	std::cout<<"Ground truth  "<<" -> "<<std::setw(15)<<std::setprecision(8)
			<<std::scientific<<groundtruth(rmin,rmax)<<std::endl;
	std::cout<<std::endl<<"Basic"<<std::endl;	
	test("Trapezoidal   ",integrator_quadrature(trapezoidal), f,rmin,rmax);
	test("Simpson       ",integrator_quadrature(simpson),f,rmin,rmax);
	test("Boole         ",integrator_quadrature(boole), f,rmin,rmax);
	std::cout<<std::endl<<"Steps"<<std::endl;
	test("Trapezoidal x2",integrator_quadrature(steps<2>(trapezoidal)), f,rmin,rmax);
	test("Trapezoidal x4",integrator_quadrature(steps<4>(trapezoidal)), f,rmin,rmax);
	test("Simpson x2    ",integrator_quadrature(steps<2>(simpson)), f,rmin,rmax);
	test("Trapezoidal x8",integrator_quadrature(steps<8>(trapezoidal)), f,rmin,rmax);
	test("Simpson x4    ",integrator_quadrature(steps<4>(simpson)), f,rmin,rmax);
	test("Boole x2      ",integrator_quadrature(steps<2>(boole)), f,rmin,rmax);
	test("Trapezoid. x16",integrator_quadrature(steps<16>(trapezoidal)), f,rmin,rmax);
	test("Simpson x8    ",integrator_quadrature(steps<8>(simpson)), f,rmin,rmax);
	test("Boole x4      ",integrator_quadrature(steps<4>(boole)), f,rmin,rmax);
	test("Trapezoid. x32",integrator_quadrature(steps<32>(trapezoidal)), f,rmin,rmax);
	test("Simpson x16   ",integrator_quadrature(steps<16>(simpson)), f,rmin,rmax);
	test("Boole x8      ",integrator_quadrature(steps<8>(boole)), f,rmin,rmax);


	std::cout<<std::endl<<"Adaptive tolerance"<<std::endl;
	test("Tr.Simp. e-2  ",integrator_adaptive_tolerance(nested(simpson,trapezoidal),1.e-2),f,rmin,rmax);
	test("Tr.Simp. e-4  ",integrator_adaptive_tolerance(nested(simpson,trapezoidal),1.e-4),f,rmin,rmax);
	test("Tr.Simp. e-6  ",integrator_adaptive_tolerance(nested(simpson,trapezoidal),1.e-6),f,rmin,rmax);
	test("Simp.Boole e-2",integrator_adaptive_tolerance(nested(boole,simpson),1.e-2),f,rmin,rmax);
	test("Simp.Boole e-4",integrator_adaptive_tolerance(nested(boole,simpson),1.e-4),f,rmin,rmax);
	test("Simp.Boole e-6",integrator_adaptive_tolerance(nested(boole,simpson),1.e-6),f,rmin,rmax);

	std::cout<<std::endl<<"Adaptive iterations"<<std::endl;
	test("Tr.Simp.    1 ",integrator_adaptive_iterations(nested(simpson,trapezoidal),1),f,rmin,rmax);
	test("Tr.Simp.   10 ",integrator_adaptive_iterations(nested(simpson,trapezoidal),10),f,rmin,rmax);
	test("Tr.Simp.  100 ",integrator_adaptive_iterations(nested(simpson,trapezoidal),100),f,rmin,rmax);
	test("Tr.Simp. 1000 ",integrator_adaptive_iterations(nested(simpson,trapezoidal),1000),f,rmin,rmax);
	test("Simp.Boole   1",integrator_adaptive_iterations(nested(boole,simpson),1),f,rmin,rmax);
	test("Simp.Boole  10",integrator_adaptive_iterations(nested(boole,simpson),10),f,rmin,rmax);
	test("Simp.Boole 100",integrator_adaptive_iterations(nested(boole,simpson),100),f,rmin,rmax);
	test("Simp.Boole1000",integrator_adaptive_iterations(nested(boole,simpson),1000),f,rmin,rmax);


	std::cout<<std::endl<<"Monte Carlo"<<std::endl;
	test("Uniform     1 ",integrator_monte_carlo_uniform(1),f,rmin,rmax);
	test("Uniform    10 ",integrator_monte_carlo_uniform(10),f,rmin,rmax);
	test("Uniform   100 ",integrator_monte_carlo_uniform(100),f,rmin,rmax);
	test("Uniform  1000 ",integrator_monte_carlo_uniform(1000),f,rmin,rmax);

	return 0;
}



