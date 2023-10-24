#include "../viltrum.h"
#include <iostream>
#include <iomanip>


using namespace viltrum;

class Function {
public:
    template<std::size_t N>
    float operator()(const std::array<float,N>& x) const {
        float r = 1;
        for (auto xi : x) r*=xi;
        return r;
    }
};


template<typename R>
void test_region_integral(const char* name, const R& r, const std::array<float,4>& range_min, const std::array<float,4>& range_max) {
    std::cout<<name<<"\t"<<"  Integral  \t      SUB1\t      SUB2\t      SUB3\t       ALL"<<std::endl;
    std::cout<<"Default   -  \t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral()<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange(resize<1>(range_min),resize<1>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange(resize<2>(range_min),resize<2>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange(resize<3>(range_min),resize<3>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange(range_min,range_max)<<std::endl;
    std::cout<<"First     -  \t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral()<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange_first(resize<1>(range_min),resize<1>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange_first(resize<2>(range_min),resize<2>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange_first(resize<3>(range_min),resize<3>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange_first(range_min,range_max)<<std::endl;
    std::cout<<"Last      -  \t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral()<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange_last(resize<1>(range_min),resize<1>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange_last(resize<2>(range_min),resize<2>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange_last(resize<3>(range_min),resize<3>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<r.integral_subrange_last(range_min,range_max)<<std::endl;
	auto p = r.polynomial();
	std::cout<<"Polynomial-  \t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<p.integral()<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<p.integral(resize<1>(range_min),resize<1>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<p.integral(resize<2>(range_min),resize<2>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<p.integral(resize<3>(range_min),resize<3>(range_max))<<"\t";
    std::cout<<std::fixed<<std::setprecision(3)<<std::setw(10)<<p.integral(range_min,range_max)<<std::endl;


    std::cout<<std::endl;
}

template<typename R>
void test_region_polynomial(const char* name, const R& r, const std::array<float,4>& range_min, const std::array<float,4>& range_max) {
	
	auto p = r.polynomial();
	std::cout<<name<<" - polynomial approximation"<<std::endl;
	std::cout<<"Region values : "<<r.approximation_at(range_min)<<"\t"<<r.approximation_at(range_max)<<"\t"<<r.approximation_at(std::array{0.0f,0.0f,0.0f,0.0f})<<"\t"<<r.approximation_at(std::array{0.5f,-0.5f,0.5f,-0.5f})<<std::endl;
	std::cout<<"Polyn. values : "<<p(range_min)<<"\t"<<p(range_max)<<"\t"<<p(std::array{0.0f,0.0f,0.0f,0.0f})<<"\t"<<p(std::array{0.5f,-0.5f,0.5f,-0.5f})<<std::endl;
}

int main(int argc, char **argv) {
    Function f;
    std::array<float,4> range_min{-1,-1,-1,-1};
    std::array<float,4> range_max{ 2, 2, 2, 2};
	
    auto region_trapezoidal = region(f,trapezoidal,range_min,range_max);   
    auto region_simpson     = region(f,simpson,range_min,range_max);   
    auto region_boole       = region(f,boole,range_min,range_max);
	
    test_region_polynomial("Trapezoidal",region_trapezoidal,range_min,range_max);
    test_region_polynomial("Simpson    ",region_simpson,range_min,range_max);
    test_region_polynomial("Boole      ",region_boole,range_min,range_max);
	
    test_region_integral("Trapezoidal",region_trapezoidal,range_min,range_max);
    test_region_integral("Simpson    ",region_simpson,range_min,range_max);
    test_region_integral("Boole      ",region_boole,range_min,range_max);
}

