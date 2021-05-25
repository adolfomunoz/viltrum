#include "../viltrum.h"
#include <iostream>
#include <iomanip>
#include <functional>
#include <numeric>

using namespace viltrum;

class Function {
	mutable unsigned long evals = 0;
public:
    template<typename Float, std::size_t DIM>
	auto operator()(const std::array<Float, DIM>& pos) const {
		++evals;
		return std::accumulate(pos.begin(),pos.end(),Float(0));
	}
	
	unsigned long evaluations() const { return evals; }
};

template<std::size_t DIM>
Range<double,DIM> default_range() {
	std::array<double,DIM> a, b;
	std::fill(a.begin(),a.end(),0.0f);
	std::fill(b.begin(),b.end(),1.0f);
	return Range<double,DIM>(a,b);
}

template<std::size_t DIM, typename Stepper>
long samples_per_iteration(const Stepper& stepper) {
	Function fant, f;
	/*auto sant = */integrator_stepper(stepper,1).integrate(fant,default_range<DIM>());
	/*auto s = */integrator_stepper(stepper,2).integrate(f,default_range<DIM>());
	return f.evaluations()-fant.evaluations();
}

template<typename Stepper>
void test(const char* name, const Stepper& stepper) {
	std::cout<<name<<"\t";
	std::cout<<std::setw(10)<<samples_per_iteration<1>(stepper);
	std::cout<<std::setw(10)<<samples_per_iteration<2>(stepper);
	std::cout<<std::setw(10)<<samples_per_iteration<3>(stepper);
	std::cout<<std::setw(10)<<samples_per_iteration<4>(stepper);
	std::cout<<std::setw(10)<<samples_per_iteration<5>(stepper);
	std::cout<<std::setw(10)<<samples_per_iteration<6>(stepper);
	std::cout<<std::endl;
}

int main(int argc, char **argv) {
	std::cout<<"   DIMENSIONS\t";
	for (int i = 1; i<=6; ++i ) std::cout<<std::setw(10)<<i;
	std::cout<<std::endl;
	test("Monte-Carlo  ",stepper_monte_carlo_uniform());
	test("Simpson-Trap.",stepper_adaptive(nested(simpson,trapezoidal)));
	test("Boole-Simpson",stepper_adaptive(nested(boole,simpson)));
}
